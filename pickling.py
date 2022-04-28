# pickling.py: Pickle and Unpickle entire NRPy+ environment.
#              This is a key element to NRPy+ parallel codegen.

# Author: Zachariah B. Etienne
#         zachetie **at** gmail **dot* com

from outputC import outC_function_dict, outC_function_prototype_dict, outC_function_outdir_dict  # NRPy+: Core C code output module
import NRPy_param_funcs as par   # NRPy+: Parameter interface
import grid as gri               # NRPy+: Functions having to do with numerical grids

def pickle_NRPy_env():
    # Store all NRPy+ environment variables to an output string so NRPy+ environment from within this subprocess can be easily restored
    import pickle
    # https://www.pythonforthelab.com/blog/storing-binary-data-and-serializing/
    outstr = []
    outstr.append(pickle.dumps(len(gri.glb_gridfcs_list)))
    for lst in gri.glb_gridfcs_list:
        outstr.append(pickle.dumps(lst.gftype))
        outstr.append(pickle.dumps(lst.name))
        outstr.append(pickle.dumps(lst.rank))
        outstr.append(pickle.dumps(lst.DIM))

    outstr.append(pickle.dumps(len(par.glb_params_list)))
    for lst in par.glb_params_list:
        outstr.append(pickle.dumps(lst.type))
        outstr.append(pickle.dumps(lst.module))
        outstr.append(pickle.dumps(lst.parname))
        outstr.append(pickle.dumps(lst.defaultval))

    outstr.append(pickle.dumps(len(par.glb_Cparams_list)))
    for lst in par.glb_Cparams_list:
        outstr.append(pickle.dumps(lst.type))
        outstr.append(pickle.dumps(lst.module))
        outstr.append(pickle.dumps(lst.parname))
        outstr.append(pickle.dumps(lst.defaultval))

    outstr.append(pickle.dumps(len(outC_function_dict)))
    for Cfuncname, Cfunc in outC_function_dict.items():
        outstr.append(pickle.dumps(Cfuncname))
        outstr.append(pickle.dumps(Cfunc))

    outstr.append(pickle.dumps(len(outC_function_prototype_dict)))
    for Cfuncname, Cfuncprototype in outC_function_prototype_dict.items():
        outstr.append(pickle.dumps(Cfuncname))
        outstr.append(pickle.dumps(Cfuncprototype))

    outstr.append(pickle.dumps(len(outC_function_outdir_dict)))
    for Cfuncname, Cfuncoutdir in outC_function_outdir_dict.items():
        outstr.append(pickle.dumps(Cfuncname))
        outstr.append(pickle.dumps(Cfuncoutdir))

    return outstr

def unpickle_NRPy_env(NRPyEnvVars):
    import pickle
    # https://www.pythonforthelab.com/blog/storing-binary-data-and-serializing/
    grfcs_list = []
    param_list = []
    Cparm_list = []

    outCfunc_dict = {}
    outCfuncproto_dict = {}
    outCfuncoutdir_dict = {}

    for WhichParamSet in NRPyEnvVars[0]:
        # gridfunctions
        i=0
        num_elements = pickle.loads(WhichParamSet[i]); i+=1
        for lst in range(num_elements):
            grfcs_list.append(gri.glb_gridfc(gftype=pickle.loads(WhichParamSet[i+0]),
                                         name  =pickle.loads(WhichParamSet[i+1]),
                                         rank  =pickle.loads(WhichParamSet[i+2]),
                                         DIM   =pickle.loads(WhichParamSet[i+3]))) ; i+=4
        # parameters
        num_elements = pickle.loads(WhichParamSet[i]); i+=1
        for lst in range(num_elements):
            param_list.append(par.glb_param(type      =pickle.loads(WhichParamSet[i+0]),
                                            module    =pickle.loads(WhichParamSet[i+1]),
                                            parname   =pickle.loads(WhichParamSet[i+2]),
                                            defaultval=pickle.loads(WhichParamSet[i+3]))) ; i+=4
        # Cparameters
        num_elements = pickle.loads(WhichParamSet[i]); i+=1
        for lst in range(num_elements):
            Cparm_list.append(par.glb_Cparam(type      =pickle.loads(WhichParamSet[i+0]),
                                             module    =pickle.loads(WhichParamSet[i+1]),
                                             parname   =pickle.loads(WhichParamSet[i+2]),
                                             defaultval=pickle.loads(WhichParamSet[i+3]))) ; i+=4
        # outC_func_dict
        num_elements = pickle.loads(WhichParamSet[i]); i+=1
        for lst in range(num_elements):
            funcname = pickle.loads(WhichParamSet[i+0])
            funcbody = pickle.loads(WhichParamSet[i+1]) ; i+=2
            outCfunc_dict[funcname] = funcbody

        # outC_function_prototype_dict
        num_elements = pickle.loads(WhichParamSet[i]); i+=1
        for lst in range(num_elements):
            funcname  = pickle.loads(WhichParamSet[i+0])
            funcproto = pickle.loads(WhichParamSet[i+1]) ; i+=2
            outCfuncproto_dict[funcname] = funcproto

        # outC_function_outdir_dict
        num_elements = pickle.loads(WhichParamSet[i]); i+=1
        for lst in range(num_elements):
            funcname   = pickle.loads(WhichParamSet[i+0])
            funcoutdir = pickle.loads(WhichParamSet[i+1]) ; i+=2
            outCfuncoutdir_dict[funcname] = funcoutdir

    grfcs_list_uniq = []
    for gf_ntuple_stored in grfcs_list:
        found_gf = False
        for gf_ntuple_new in grfcs_list_uniq:
            if gf_ntuple_new == gf_ntuple_stored:
                found_gf = True
        if found_gf == False:
            grfcs_list_uniq.append(gf_ntuple_stored)

    param_list_uniq = []
    for pr_ntuple_stored in param_list:
        found_pr = False
        for pr_ntuple_new in param_list_uniq:
            if pr_ntuple_new == pr_ntuple_stored:
                found_pr = True
        if found_pr == False:
            param_list_uniq.append(pr_ntuple_stored)

    # Set glb_paramsvals_list:
    # Step 1: Reset all paramsvals to their defaults
    #  BAD IDEA: OVERWRITTEN DEFAULTS SHOULD BE KEPT.
    # par.glb_paramsvals_list = []
    # for parm in param_list_uniq:
    #     par.glb_paramsvals_list.append(parm.defaultval)

    Cparm_list_uniq = []
    for Cp_ntuple_stored in Cparm_list:
        found_Cp = False
        for Cp_ntuple_new in Cparm_list_uniq:
            if Cp_ntuple_new == Cp_ntuple_stored:
                found_Cp = True
        if found_Cp == False:
            Cparm_list_uniq.append(Cp_ntuple_stored)

    gri.glb_gridfcs_list = []
    par.glb_params_list  = []
    par.glb_Cparams_list = []

    gri.glb_gridfcs_list = grfcs_list_uniq
    par.glb_params_list  = param_list_uniq
    par.glb_Cparams_list = Cparm_list_uniq
    for key, item in outCfunc_dict.items():
        outC_function_dict[key] = item

    for key, item in outCfuncproto_dict.items():
        outC_function_prototype_dict[key] = item

    for key, item in outCfuncoutdir_dict.items():
        outC_function_outdir_dict[key] = item
