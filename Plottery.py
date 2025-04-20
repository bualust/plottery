#! /usr/bin/python3
import sys
import os
import glob
import uproot
import boost_histogram as bh
import matplotlib.pyplot as plt
import awkward as ak
import mplhep as hep
import numpy as np
import yaml
from matplotlib.gridspec import GridSpec
from optparse import OptionParser
import re

# ==============================================================
#  SETTINGS
# ==============================================================
hep.style.use(hep.style.ATLAS)
Variables = []
VariablesTitle = []
Varriables_nbins = []
VariablesBinLowEdge = []
VariablesBinHighEdge = []
VariablesBinning = []

# =================================================================================
#  main
# =================================================================================
def main():
    parser = OptionParser()
    parser.add_option("--config", dest="config_file", default="config.yaml",
                            help="Pass config file [default: %default]", metavar="STRING")

    try:
        (options, args) = parser.parse_args()
    except:
        parser.print_help()
        exit()

    config_file = options.config_file
    info('Opening config file: '+config_file)
    with open(config_file) as yaml_file:
        yaml_data = yaml.load(yaml_file, Loader=yaml.FullLoader)

    plotting(yaml_data)
    exit()

# ==============================================================
#  Opens ROOT files with uproot and gets the branches/Variables
# ==============================================================
def getdata(yaml_data,mc_weight,selection_cuts,treename):

    allgood("Preparing to read input files")
    #If casting is used more files are considered for
    #a single process so they are grouped in FilesToOpen
    FilesToOpen = []
    for proc in yaml_data['Processes']['File']:
        filenames = glob.glob(yaml_data['InputDirectory']+proc+'*')
        append_tree = []
        for ff in filenames: append_tree.append(ff+':'+treename)
        FilesToOpen.append(append_tree)

    #convernting yaml's Variables into arrays
    get_variables(yaml_data['Variables'])

    #adding mc_weight
    if mc_weight!='':
        Variables.append(mc_weight)
        VariablesTitle.append(mc_weight)

    #adding variables needed for selections
    variables = Variables
    variables = list_of_variables(variables,selection_cuts,yaml_data)

    # deal with selection on vectors
    sel_var_vectors = get_sel_vectors(variables)

    #data will be in the form of a dictionary
    #with a key assigned to each process
    data = {}
    for (myfiles,proc) in zip(FilesToOpen,yaml_data['Processes']['Name']):
        if proc=='Data': variables = [item for item in variables if item != mc_weight]
        var_arr = uproot.concatenate(myfiles, variables, allow_missing=True)
        variables_arr   = {var:var_arr[var] for var in variables}
        data[proc] = variables_arr
        if proc=='Data': variables.append(mc_weight)

    #making vectors flat
    #first from variables to be plotted
    var_to_remove = []
    for proc in yaml_data['Processes']['Name']:
        process_dic = data[proc]
        yaml_var = yaml_data['Variables']
        for i in range(len(yaml_var)):
            if yaml_var[i]['Index']>=0:
                flat_vector(process_dic, yaml_var[i]['Index'], yaml_var[i]['branch_name'])
                var_to_remove.append(yaml_var[i]['branch_name'])
        data[proc] = process_dic

    #then from selection variables
    for proc in yaml_data['Processes']['Name']:
        process_dic = data[proc]
        sel_vars = sel_var_vectors.keys()
        for i_sel_vars in sel_vars:
            flat_vector(process_dic, sel_var_vectors[i_sel_vars] , i_sel_vars)
            var_to_remove.append(i_sel_vars)
        data[proc] = process_dic

    #and removing vectors of flattened variables
    var_to_remove = np.unique(var_to_remove)
    for proc in yaml_data['Processes']['Name']:
        for var_rm in var_to_remove: del data[proc][var_rm]

    data = select(data,yaml_data,selection_cuts,sel_var_vectors)
    allgood('Dictionary with data is ready!')
    return data

# ==============================================================
# takes list of variables to retrieve from root file
# if one of the is a vector from a selection it strips [*]
# from its name and return dictiornary with original index
# ==============================================================
def get_sel_vectors(variables):
    sel_var_vectors = {}
    var_index = 0
    for var in variables:
        if "[" in var:
            vec_var_ind = re.findall(r"\[\s*\+?(-?\d+)\s*\]", var)
            int_vec_var_ind = int(vec_var_ind[0])
            var = re.sub(r"\[\s*\+?(-?\d+)\s*\]", "", var)
            var = var.replace("[","")
            var = var.replace("]","")
            variables[var_index] = var
            sel_var_vectors[var] = int_vec_var_ind
        var_index = var_index+1
    return sel_var_vectors

# ==============================================================
#  adds to the list of variables to get with uproot from selections
# ==============================================================
def list_of_variables(variables,selection_cuts,yaml_data):

    #first loading variables from sample selection
    for sample_sel in yaml_data["Processes"]["Selection"]:
        if sample_sel=="": continue
        new_list = selection_to_variables(sample_sel)
        for s in new_list: variables.append(s)

    #then loading variables from events selection
    if selection_cuts == "": return variables
    else:
        new_list = selection_to_variables(selection_cuts)
        for s in new_list: variables.append(s)


    return variables

# ==============================================================
#  strips the selection, operators and cut values from the
#  selections strings and builds a list of variables to be loaded
#  from the input files
# ==============================================================
def selection_to_variables(selection):
    sample_sel = re.sub("abs|not","",selection)
    sample_sel = re.sub(r"\(|\)","",sample_sel)
    pattern = r"[\s><=!]+(?:and|or)?[\s><=!0-9]*\d*(?:\.\d+)?"
    new_string = re.sub(pattern, ',',sample_sel)  # replace matched substrings with comma
    new_string = re.sub(r'(and|or),', ',', new_string)  # remove any remaining and/or followed by comma
    new_list = new_string.split(',')  # split the string at each comma
    new_list = [s.strip() for s in new_list]  # remove leading/trailing whitespaces
    new_list = [s for s in new_list if s != '' and s not in ['and', 'or']]  # remove empty strings, and 'and'/'or' from the list
    return new_list

# ==============================================================
#  flatten vectors variables
# ==============================================================
def flat_vector(process_dic, index, name):
    var_arr = []
    for ll in process_dic[name]: var_arr.append(ll[index])
    var_arr = ak.Array(var_arr)
    new_name = name+"_"+str(index)
    process_dic[new_name] = var_arr
    for i in range(len(Variables)):
        if Variables[i]==name: Variables[i]=new_name; break

# ==============================================================
#  Convert yaml variables into arrays
# ==============================================================
def get_variables(yaml_var):

    for i in range(len(yaml_var)):
        Variables.append(yaml_var[i]['branch_name'])
        VariablesTitle.append(yaml_var[i]['Title'])
        Varriables_nbins.append(yaml_var[i]['Nbins'])
        VariablesBinLowEdge.append(yaml_var[i]['LowerRange'])
        VariablesBinHighEdge.append(yaml_var[i]['UpperRange'])
        VariablesBinning.append(yaml_var[i]['Binning'])

    return allgood('Variables yaml to array conversion: DONE')

# ==============================================================
#  Apply event selection
# ==============================================================
def select(data,yaml_data,selection_cuts,sel_var_vectors):

    info('Applying event selection')
    for proc,sel in zip(yaml_data['Processes']['Name'],yaml_data['Processes']['Selection']):
        process_dic = data[proc]
        #Switching to pandas to apply selection
        pandas_proc_dic = ak.to_pandas(process_dic)
        sample_selection_cuts=selection_cuts+sel
        if "[" in sample_selection_cuts:
            sample_selection_cuts = sample_selection_cuts.replace("[","_")
            sample_selection_cuts = sample_selection_cuts.replace("]","")
        after_sel_dic = pandas_proc_dic.query(sample_selection_cuts)
        data[proc] = after_sel_dic

    return data

# ==============================================================
#  Building histograms, orvelay them in same figure and save
# ==============================================================
def plotting(yaml_data):

    mc_weight      = yaml_data['mc_weight']
    OutDir         = yaml_data['OutDir']
    OutDir         = 'Output/'+OutDir+'/'
    label          = yaml_data['Label']
    lumi           = yaml_data['Lumi']
    normalise      = yaml_data['Normalise']
    plot_label     = '$\sqrt{s}=13$ TeV, '+str(lumi)+' fb$^{-1}$\n'+label
    selection_cuts = yaml_data['SelectionCuts']
    treename       = yaml_data['TreeName']

    data = getdata(yaml_data,mc_weight,selection_cuts,treename)
    allgood('Plotting time!')

    for var,var_name,var_nbins,var_low,var_high,VarBinning in zip(Variables,VariablesTitle,Varriables_nbins,VariablesBinLowEdge,VariablesBinHighEdge,VariablesBinning):
         #if mc_weight in var: continue ## not really interested in weight_nominal histo
         fig = plt.figure(figsize=(7.0,7.0),dpi=400)
         gs = GridSpec(2,1, height_ratios=[6,1],hspace=0.1)
         main = fig.add_subplot(gs[0])
         ratio= fig.add_subplot(gs[1],sharex=main)
         histos = {}
         proc_yields = []
         for proc,rescale in zip(yaml_data['Processes']['Name'],yaml_data['Processes']['Rescale']):
            h_ax = get_xaxis(VarBinning,var_nbins,var_low,var_high)
            hist = bh.Histogram(h_ax,storage=bh.storage.Weight())
            info('Getting dictionary for process:   '+proc)
            process_dic = data[proc]
            info('Getting '+var+' values for '+proc)
            process_var = process_dic[var]
            if mc_weight=='' or proc=='Data': process_weight=1.
            else: process_weight = process_dic[mc_weight]
            lumi_scale = lumi
            if proc=='Data': lumi_scale=1
            lumi_scale=lumi_scale*rescale
            hist.fill(process_var, weight=process_weight*lumi_scale)
            this_yield = 0
            if mc_weight=='' or proc=='Data': 
                for ev in process_var: this_yield=this_yield+lumi_scale
            else:
                for weight in process_weight: this_yield=this_yield+weight*lumi_scale
            if proc=='Data': 
                hep.histplot(hist,yerr=True,histtype='errorbar',ax=main,density=normalise,color="black",marker="o",markersize=5,label=proc+" ("+str(int(this_yield))+")")
            else: 
                hep.histplot(hist,yerr=True,histtype='step',ax=main,density=normalise,label=proc+" ("+str(int(this_yield))+")")
            hep.atlas.label(data=True, lumi=lumi, label='Internal', rlabel=plot_label,ax=main)
            histos[proc] = hist
            proc_yields.append(this_yield)
         ii = 0
         for i_hist in histos:
            num_hist = histos[i_hist]
            ref_sample = yaml_data['Processes']['Name'][0]
            ref_hist = histos[ref_sample]
            if len(yaml_data['Processes']['Name'])==1: vals_ratio=num_hist.values()
            else:
                if normalise: vals_ratio = (num_hist.values()/ref_hist.values())*(proc_yields[0]/proc_yields[ii])
                else: vals_ratio = (num_hist.values()/ref_hist.values())
            if proc=='Data': 
                hep.histplot(vals_ratio, bins=ref_hist.axes[0].edges, yerr=False,histtype='step',ax=ratio,color="black")
            else: 
                hep.histplot(vals_ratio, bins=ref_hist.axes[0].edges, yerr=False,histtype='step',ax=ratio)
            ii = ii +1
         plt.xlabel(var_name)
         ymin, ymax = main.axes.get_ylim()
         main.set_ylim(ymin,ymax+(ymax-ymin)*0.3)
         ratio.set_ylim(0.0,2.0)
         main.axes.xaxis.set_visible(False)
         main.legend(loc='upper right')
         info('Saving plots in '+OutDir)
         os.system('mkdir -p '+OutDir)
         ratio.axhline(y=1.0, color='grey',linestyle='--')
         fig.savefig(OutDir+var+'.pdf')

    return allgood('All plots saved')

# ===========================================================
#  Builds histogram xaxis binning
# ===========================================================
def get_xaxis(VarBinning,var_nbins,var_low,var_high):

    if(VarBinning==-99):
        h_ax = bh.axis.Regular(var_nbins,var_low,var_high)
    else:
        floats_list = [float(item) for item in VarBinning.split(',')]
        h_ax = bh.axis.Variable(floats_list)

    return h_ax

# ===========================================================
#  Colours to throw scary messages on screen
# ===========================================================
class bcolors:
    INFO = '\033[95m' #PURPLE
    OK = '\033[92m' #GREEN
    WARNING = '\033[93m' #YELLOW
    FAIL = '\033[91m' #RED
    RESET = '\033[0m' #RESET COLOR

def warn(message): print(bcolors.WARNING+"=== "+message+" ==="+bcolors.RESET)
def fail(message): print(bcolors.FAIL+"=== "+message+" ==="+bcolors.RESET)
def allgood(message): print(bcolors.OK+"=== "+message+" ==="+bcolors.RESET)
def info(message): print(bcolors.INFO+"=== "+message+" ==="+bcolors.RESET)

# ===========================================================
# __main__
# ===========================================================
if __name__ == '__main__':
  main()
