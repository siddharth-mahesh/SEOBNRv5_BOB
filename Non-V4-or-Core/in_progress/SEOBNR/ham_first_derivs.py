# NRPy+ code to generate first derivatives of the SEOBNRv3 Hamiltonian from a list of numerical expressions computing
# said Hamiltonian. Originally written by Zach Etienne; edited and commented by Tyler Knowles.

from outputC import outputC,lhrh,superfast_uniq  # NRPy+: Core C code output module
import sympy as sp      # SymPy: The Python computer algebra package upon which NRPy+ depends
import sys              # Python module for multiplatform OS-related functions

# As of April 2021, "sp.sympify("Q+1")" fails because Q is a reserved keyword.
#   This is the workaround, courtesy Ken Sible.
custom_global_dict = {}
exec('from sympy import *', custom_global_dict)
del custom_global_dict['Q']
custom_parse_expr = lambda expr: sp.parse_expr(expr, global_dict=custom_global_dict)

# simplify_deriv() simplifies derivative expressions by removing terms equal to zero.
def simplify_deriv(lhss_deriv, rhss_deriv):
    # Create 'simp' arrays to store and manipulate derivative expressions.
    lhss_deriv_simp = []
    rhss_deriv_simp = []
    # Append terms to 'simp' arrays.
    for i in range(len(rhss_deriv)):
        lhss_deriv_simp.append(lhss_deriv[i])
        rhss_deriv_simp.append(rhss_deriv[i])
    # For each term equal to zero, loop through all expressions and replace that variable with the number zero.
    for i in range(len(rhss_deriv_simp)):
        if rhss_deriv_simp[i] == 0:
            for j in range(i + 1, len(rhss_deriv_simp)):
                for var in rhss_deriv_simp[j].free_symbols:
                    if str(var) == str(lhss_deriv_simp[i]):
                        rhss_deriv_simp[j] = rhss_deriv_simp[j].subs(var, 0)
    # Create 'zero' array to store terms to be removed from derivative expressions.
    zero_elements_to_remove = []
    # Loop over all terms and add those equal to zero to 'zero' array.
    for i in range(len(rhss_deriv_simp)):
        if rhss_deriv_simp[i] == sp.sympify(0):
            zero_elements_to_remove.append(i)
    count = 0
    # Remove from derivative list all elements of 'zero' array.
    for i in range(len(zero_elements_to_remove)):
        del lhss_deriv_simp[zero_elements_to_remove[i] + count]
        del rhss_deriv_simp[zero_elements_to_remove[i] + count]
        count -= 1
    # Return simplified derivative expressions.
    return lhss_deriv_simp, rhss_deriv_simp


# deriv_onevar() replaces variable derivatives with 1 or 0 depending on which partial derivaitve is computed.  For
# example, pass 'xprm=1' to replace each instance of 'xprm' with 1 and 'qprm' with 0 for each q in (y,z,p1,p2,p3,S1x,
# S1y,S1z,S2x,S2y,S2z).  This produces expressions which compute the partial derivative of the Hamiltonian with respect
# to x.
def deriv_onevar(lhss_deriv, rhss_deriv, xprm=0, yprm=0, zprm=0, p1prm=0, p2prm=0, p3prm=0, S1xprm=0, S1yprm=0,
                 S1zprm=0, S2xprm=0, S2yprm=0, S2zprm=0):
    if xprm + yprm + zprm + p1prm + p2prm + p3prm + S1xprm + S1yprm + S1zprm + S2xprm + S2yprm + S2zprm != 1:
        print("deriv_onevar() cannot take more than one derivative at a time!")
        sys.exit()

    # Create 'new' arrays to store and manipulate derivative terms.
    lhss_deriv_new = []
    rhss_deriv_new = []
    # Append derivative terms to 'new' arrays
    for i in range(len(rhss_deriv)):
        lhss_deriv_new.append(lhss_deriv[i])
        rhss_deriv_new.append(rhss_deriv[i])
    # Replace each instance of 'qprm', q in (x,y,z,p1,p2,p3,S1x,S1y,S1z,S2x,S2y,S2z), with either 0 or 1.
    for i in range(len(rhss_deriv_new)):
        for var in rhss_deriv_new[i].free_symbols:
            if str(var) == "xprm":
                rhss_deriv_new[i] = rhss_deriv_new[i].subs(var, xprm)
            elif str(var) == "yprm":
                rhss_deriv_new[i] = rhss_deriv_new[i].subs(var, yprm)
            elif str(var) == "zprm":
                rhss_deriv_new[i] = rhss_deriv_new[i].subs(var, zprm)
            elif str(var) == "p1prm":
                rhss_deriv_new[i] = rhss_deriv_new[i].subs(var, p1prm)
            elif str(var) == "p2prm":
                rhss_deriv_new[i] = rhss_deriv_new[i].subs(var, p2prm)
            elif str(var) == "p3prm":
                rhss_deriv_new[i] = rhss_deriv_new[i].subs(var, p3prm)
            elif str(var) == "S1xprm":
                rhss_deriv_new[i] = rhss_deriv_new[i].subs(var, S1xprm)
            elif str(var) == "S1yprm":
                rhss_deriv_new[i] = rhss_deriv_new[i].subs(var, S1yprm)
            elif str(var) == "S1zprm":
                rhss_deriv_new[i] = rhss_deriv_new[i].subs(var, S1zprm)
            elif str(var) == "S2xprm":
                rhss_deriv_new[i] = rhss_deriv_new[i].subs(var, S2xprm)
            elif str(var) == "S2yprm":
                rhss_deriv_new[i] = rhss_deriv_new[i].subs(var, S2yprm)
            elif str(var) == "S2zprm":
                rhss_deriv_new[i] = rhss_deriv_new[i].subs(var, S2zprm)
    # Simplify the derivative expressions with simplify_deriv().
    lhss_deriv_simp, rhss_deriv_simp = simplify_deriv(lhss_deriv_new, rhss_deriv_new)
    # Return simplified derivative expression.
    return lhss_deriv_simp, rhss_deriv_simp


# replace_numpy_funcs() replaces specific SymPy function names with the corresponding NumPy function names.
def replace_numpy_funcs(expression):
    return str(expression).replace("sqrt(", "sp.sqrt(").replace("Abs(", "sp.Abs(").replace("log(",
                                                                "sp.log(").replace("sign(", "sp.sign(")


# output_H_and_derivs() is the main wrapper function for computing the SEONBRv3 Hamiltonian H and the twelve first
# partial derivatives of H with respect to x, y, z, p1, p2, p3, S1x, S1y, S1z, S2x, S2y, S2z.
# TylerK: for now, only outputs dHdx, dHdpy, and dHdpz for initial condition root-finding!
def output_H_and_derivs():
    # Open and read the file of numerical expressions (written in SymPy syntax) computing the SEOBNRv3 Hamiltonian.
    #f = open("SEOBNR/Hamstring.txt", 'r')
    f = open("SEOBNR/SymPy_Hreal_on_bottom.txt", 'r')
    Hamstring = str(f.read())
    f.close()

    # Split Hamstring by carriage returns.
    Hamterms = Hamstring.splitlines()

    # Create 'lr' array to store each left-hand side and right-hand side of Hamstring as strings.
    lr = []
    # Loop over each line in Hamstring to separate the left- and right-hand sides.
    for i in range(len(Hamterms)):
        # Ignore lines with 2 or fewer characters and those starting with #
        if len(Hamterms[i]) > 2 and Hamterms[i][0] != "#":
            # Split each line by its equals sign.
            splitHamterms = Hamterms[i].split("=")
            # Append terms to the 'lr' array, removing spaces, "sp." prefixes, and replacing Lambda->Lamb (Lambda is a
            # protected keyword)
            lr.append(lhrh(lhs=splitHamterms[0].replace(" ", "").replace("Lambda", "Lamb"),
                           rhs=splitHamterms[1].replace(" ", "").replace("sp.", "").replace("Lambda", "Lamb")))
    # Declare the symbol 'xx', which we use to denote each left-hand side as a function
    xx = sp.Symbol('xx')
    # Create arrays to store simplified left- and right-hand expressions, as well as left-hand sides designated as
    # functions.
    func = []
    lhss = []
    rhss = []
    # Affix '(xx)' to each left-hand side as a function designation; separate and simplify left- and right-hand sides
    # of the numerical expressions.
    for i in range(len(lr)):
        func.append(sp.sympify(sp.Function(lr[i].lhs)(xx)))
        #lhss.append(sp.sympify(lr[i].lhs))
        #rhss.append(sp.sympify(lr[i].rhs))
        lhss.append(custom_parse_expr(lr[i].lhs))
        rhss.append(custom_parse_expr(lr[i].rhs))
        
    # Creat array for and generate a list of all the "free symbols" in the right-hand side expressions.
    full_symbol_list_with_dups = []
    for i in range(len(lr)):
        for var in rhss[i].free_symbols:
            full_symbol_list_with_dups.append(var)

    # Remove all duplicated "free symbols" from the right-hand side expressions.
    full_symbol_list = superfast_uniq(full_symbol_list_with_dups)

    # Declare input constants.
    m1, m2, eta, KK, k0, k1, dSO, dSS = sp.symbols("m1 m2 eta KK k0 k1 dSO dSS", real=True)
    tortoise, EMgamma = sp.symbols("tortoise EMgamma", real=True)
    input_constants = [m1, m2, eta, KK, k0, k1, dSO, dSS, tortoise, EMgamma]

    # Derivatives of input constants will always be zero, so remove them from the full_symbol_list.
    for inputconst in input_constants:
        for symbol in full_symbol_list:
            if str(symbol) == str(inputconst):
                full_symbol_list.remove(symbol)

    # Add symbols to the function list and replace right-hand side terms with their function equivalent.
    full_function_list = []
    for symb in full_symbol_list:
        func = sp.sympify(sp.Function(str(symb))(xx))
        full_function_list.append(func)
        for i in range(len(rhss)):
            for var in rhss[i].free_symbols:
                if str(var) == str(symb):
                    rhss[i] = rhss[i].subs(var, func)

    # Create left- and right-hand side 'deriv' arrays
    lhss_deriv = []
    rhss_deriv = []
    # Differentiate with respect to xx, remove '(xx)', and replace xx with 'prm' notation.
    for i in range(len(rhss)):
        #lhss_deriv.append(sp.sympify(str(lhss[i]) + "prm"))
        #newrhs = sp.sympify(
        #    str(sp.diff(rhss[i], xx)).replace("(xx)", "").replace(", xx", "prm").replace("Derivative", ""))
        lhss_deriv.append(custom_parse_expr(str(lhss[i])+"prm"))
        newrhs = custom_parse_expr(str(sp.diff(rhss[i],xx)).replace("(xx)","").replace(", xx","prm").replace("Derivative",""))
        rhss_deriv.append(newrhs)
    # Simplify derivative expressions with simplify_deriv()
    lhss_deriv_simp, rhss_deriv_simp = simplify_deriv(lhss_deriv, rhss_deriv)
    lhss_deriv = lhss_deriv_simp
    rhss_deriv = rhss_deriv_simp
    # Generate partial derivatives with respect to each of the twelve input variables
    lhss_deriv_x, rhss_deriv_x = deriv_onevar(lhss_deriv, rhss_deriv, xprm=1, yprm=0, zprm=0, p1prm=0, p2prm=0, p3prm=0,
                                              S1xprm=0, S1yprm=0, S1zprm=0, S2xprm=0, S2yprm=0, S2zprm=0)
    lhss_deriv_y, rhss_deriv_y = deriv_onevar(lhss_deriv, rhss_deriv, xprm=0, yprm=1, zprm=0, p1prm=0, p2prm=0, p3prm=0,
                                               S1xprm=0, S1yprm=0, S1zprm=0, S2xprm=0, S2yprm=0, S2zprm=0)
    lhss_deriv_z, rhss_deriv_z = deriv_onevar(lhss_deriv, rhss_deriv, xprm=0, yprm=0, zprm=1, p1prm=0, p2prm=0, p3prm=0,
                                               S1xprm=0, S1yprm=0, S1zprm=0, S2xprm=0, S2yprm=0, S2zprm=0)
    lhss_deriv_p1, rhss_deriv_p1 = deriv_onevar(lhss_deriv, rhss_deriv, xprm=0, yprm=0, zprm=0, p1prm=1, p2prm=0,
                                                p3prm=0, S1xprm=0, S1yprm=0, S1zprm=0, S2xprm=0, S2yprm=0, S2zprm=0)
    lhss_deriv_p2, rhss_deriv_p2 = deriv_onevar(lhss_deriv, rhss_deriv, xprm=0, yprm=0, zprm=0, p1prm=0, p2prm=1,
                                                p3prm=0, S1xprm=0, S1yprm=0, S1zprm=0, S2xprm=0, S2yprm=0, S2zprm=0)
    lhss_deriv_p3, rhss_deriv_p3 = deriv_onevar(lhss_deriv, rhss_deriv, xprm=0, yprm=0, zprm=0, p1prm=0, p2prm=0,
                                                p3prm=1, S1xprm=0, S1yprm=0, S1zprm=0, S2xprm=0, S2yprm=0, S2zprm=0)
    lhss_deriv_S1x, rhss_deriv_S1x = deriv_onevar(lhss_deriv, rhss_deriv, xprm=0, yprm=0, zprm=0, p1prm=0, p2prm=0,
                                                  p3prm=0, S1xprm=1, S1yprm=0, S1zprm=0, S2xprm=0, S2yprm=0, S2zprm=0)
    lhss_deriv_S1y, rhss_deriv_S1y = deriv_onevar(lhss_deriv, rhss_deriv, xprm=0, yprm=0, zprm=0, p1prm=0, p2prm=0,
                                                  p3prm=0, S1xprm=0, S1yprm=1, S1zprm=0, S2xprm=0, S2yprm=0, S2zprm=0)
    lhss_deriv_S1z, rhss_deriv_S1z = deriv_onevar(lhss_deriv, rhss_deriv, xprm=0, yprm=0, zprm=0, p1prm=0, p2prm=0,
                                                  p3prm=0, S1xprm=0, S1yprm=0, S1zprm=1, S2xprm=0, S2yprm=0, S2zprm=0)
    lhss_deriv_S2x, rhss_deriv_S2x = deriv_onevar(lhss_deriv, rhss_deriv, xprm=0, yprm=0, zprm=0, p1prm=0, p2prm=0,
                                                  p3prm=0, S1xprm=0, S1yprm=0, S1zprm=0, S2xprm=1, S2yprm=0, S2zprm=0)
    lhss_deriv_S2y, rhss_deriv_S2y = deriv_onevar(lhss_deriv, rhss_deriv, xprm=0, yprm=0, zprm=0, p1prm=0, p2prm=0,
                                                  p3prm=0, S1xprm=0, S1yprm=0, S1zprm=0, S2xprm=0, S2yprm=1, S2zprm=0)
    lhss_deriv_S2z, rhss_deriv_S2z = deriv_onevar(lhss_deriv, rhss_deriv, xprm=0, yprm=0, zprm=0, p1prm=0, p2prm=0,
                                                  p3prm=0, S1xprm=0, S1yprm=0, S1zprm=0, S2xprm=0, S2yprm=0, S2zprm=1)
    
    ## functions for initial conditions and waveforms outputted separately
    
    with open("SEOBNR_Playground_Pycodes/new_dHdx.py", "w") as file:
        file.write("""from __future__ import division
import numpy as np
def new_compute_dHdx(m1, m2, EMgamma, tortoise, x, y, z, p1, p2, p3, S1x, S1y, S1z, S2x, S2y, S2z):
""")
        for i in range(len(lr) - 1):
            file.write("    " + lr[i].lhs + " = " + str(lr[i].rhs).replace("Rational(", "np.true_divide(").replace("sqrt(", "np.sqrt(").replace("log(", "np.log(").replace("sign(", "np.sign(").replace("Abs(", "np.abs(").replace("pi", "np.pi") + "\n")
        for i in range(len(lhss_deriv_x)):
            file.write("    " + str(lhss_deriv_x[i]).replace("prm", "prm_x") + " = " + replace_numpy_funcs(rhss_deriv_x[i]).replace("prm", "prm_x").replace("sp.sqrt(","np.sqrt(").replace("sp.log(","np.log(").replace("sp.sign(","np.sign(").replace("sp.Abs(", "np.abs(") + "\n")
        file.write("    return np.array([Hrealprm_x])")
                                                  
                                                  
    with open("SEOBNR_Playground_Pycodes/new_dHdp1.py", "w") as file:
        file.write("""from __future__ import division
import numpy as np
def new_compute_dHdp1(m1, m2, EMgamma, tortoise, x, y, z, p1, p2, p3, S1x, S1y, S1z, S2x, S2y, S2z):
""")
        for i in range(len(lr) - 1):
            file.write("    " + lr[i].lhs + " = " + str(lr[i].rhs).replace("Rational(", "np.true_divide(").replace("sqrt(", "np.sqrt(").replace("log(", "np.log(").replace("sign(", "np.sign(").replace("Abs(", "np.abs(").replace("pi", "np.pi") + "\n")
        for i in range(len(lhss_deriv_p1)):
            file.write("    " + str(lhss_deriv_p1[i]).replace("prm", "prm_p1") + " = " + replace_numpy_funcs(rhss_deriv_p1[i]).replace("prm", "prm_p1").replace("sp.sqrt(","np.sqrt(").replace("sp.log(","np.log(").replace("sp.sign(","np.sign(").replace("sp.Abs(", "np.abs(") + "\n")
        file.write("    return np.array([Hrealprm_p1])")        
        
    with open("SEOBNR_Playground_Pycodes/new_dHdp2.py", "w") as file:
        file.write("""from __future__ import division
import numpy as np
def new_compute_dHdp2(m1, m2, EMgamma, tortoise, x, y, z, p1, p2, p3, S1x, S1y, S1z, S2x, S2y, S2z):
""")
        for i in range(len(lr) - 1):
            file.write("    " + lr[i].lhs + " = " + str(lr[i].rhs).replace("Rational(", "np.true_divide(").replace("sqrt(", "np.sqrt(").replace("log(", "np.log(").replace("sign(", "np.sign(").replace("Abs(", "np.abs(").replace("pi", "np.pi") + "\n")
        for i in range(len(lhss_deriv_p2)):
            file.write("    " + str(lhss_deriv_p2[i]).replace("prm", "prm_p2") + " = " + replace_numpy_funcs(rhss_deriv_p2[i]).replace("prm", "prm_p2").replace("sp.sqrt(","np.sqrt(").replace("sp.log(","np.log(").replace("sp.sign(","np.sign(").replace("sp.Abs(", "np.abs(") + "\n")
        file.write("    return np.array([Hrealprm_p2])")

    with open("SEOBNR_Playground_Pycodes/new_dHdp3.py", "w") as file:
        file.write("""from __future__ import division
import numpy as np
def new_compute_dHdp3(m1, m2, EMgamma, tortoise, x, y, z, p1, p2, p3, S1x, S1y, S1z, S2x, S2y, S2z):
""")
        for i in range(len(lr) - 1):
            file.write("    " + lr[i].lhs + " = " + str(lr[i].rhs).replace("Rational(", "np.true_divide(").replace("sqrt(", "np.sqrt(").replace("log(", "np.log(").replace("sign(", "np.sign(").replace("Abs(", "np.abs(").replace("pi", "np.pi") + "\n")
        for i in range(len(lhss_deriv_p3)):
            file.write("    " + str(lhss_deriv_p3[i]).replace("prm", "prm_p3") + " = " + replace_numpy_funcs(rhss_deriv_p3[i]).replace("prm", "prm_p3").replace("sp.sqrt(","np.sqrt(").replace("sp.log(","np.log(").replace("sp.sign(","np.sign(").replace("sp.Abs(", "np.abs(") + "\n")
        file.write("    return np.array([Hrealprm_p3])")

    
    ## functions for the integrations go straight into one file
    
    with open("SEOBNR_Playground_Pycodes/new_dH.py", "w") as file:
        file.write("""from __future__ import division
import numpy as np
def new_compute_dHdx(m1, m2, EMgamma, tortoise, x, y, z, p1, p2, p3, S1x, S1y, S1z, S2x, S2y, S2z):
""")
        for i in range(len(lr) - 1):
            file.write("    " + lr[i].lhs + " = " + str(lr[i].rhs).replace("Rational(", "np.true_divide(").replace("sqrt(", "np.sqrt(").replace("log(", "np.log(").replace("sign(", "np.sign(").replace("Abs(", "np.abs(").replace("pi", "np.pi") + "\n")
        for i in range(len(lhss_deriv_x)):
            file.write("    " + str(lhss_deriv_x[i]).replace("prm", "prm_x") + " = " + replace_numpy_funcs(rhss_deriv_x[i]).replace("prm", "prm_x").replace("sp.sqrt(","np.sqrt(").replace("sp.log(","np.log(").replace("sp.sign(","np.sign(").replace("sp.Abs(", "np.abs(") + "\n")
        file.write("    return np.array([Hrealprm_x, csiprm_x])\n")
    
    with open("SEOBNR_Playground_Pycodes/new_dH.py", "a") as file:
        file.write("""def new_compute_dHdy(m1, m2, EMgamma, tortoise, x, y, z, p1, p2, p3, S1x, S1y, S1z, S2x, S2y, S2z):
""")
        for i in range(len(lr) - 1):
            file.write("    " + lr[i].lhs + " = " + str(lr[i].rhs).replace("Rational(", "np.true_divide(").replace("sqrt(", "np.sqrt(").replace("log(", "np.log(").replace("sign(", "np.sign(").replace("Abs(", "np.abs(").replace("pi", "np.pi") + "\n")
        for i in range(len(lhss_deriv_y)):
            file.write("    " + str(lhss_deriv_y[i]).replace("prm", "prm_y") + " = " + replace_numpy_funcs(rhss_deriv_y[i]).replace("prm", "prm_y").replace("sp.sqrt(","np.sqrt(").replace("sp.log(","np.log(").replace("sp.sign(","np.sign(").replace("sp.Abs(", "np.abs(") + "\n")
        file.write("    return np.array([Hrealprm_y, csiprm_y])\n")

    with open("SEOBNR_Playground_Pycodes/new_dH.py", "a") as file:
        file.write("""def new_compute_dHdz(m1, m2, EMgamma, tortoise, x, y, z, p1, p2, p3, S1x, S1y, S1z, S2x, S2y, S2z):
""")
        for i in range(len(lr) - 1):
            file.write("    " + lr[i].lhs + " = " + str(lr[i].rhs).replace("Rational(", "np.true_divide(").replace("sqrt(", "np.sqrt(").replace("log(", "np.log(").replace("sign(", "np.sign(").replace("Abs(", "np.abs(").replace("pi", "np.pi") + "\n")
        for i in range(len(lhss_deriv_z)):
            file.write("    " + str(lhss_deriv_z[i]).replace("prm", "prm_z") + " = " + replace_numpy_funcs(rhss_deriv_y[i]).replace("prm", "prm_z").replace("sp.sqrt(","np.sqrt(").replace("sp.log(","np.log(").replace("sp.sign(","np.sign(").replace("sp.Abs(", "np.abs(") + "\n")
        file.write("    return np.array([Hrealprm_z, csiprm_z])\n")
        
    with open("SEOBNR_Playground_Pycodes/new_dH.py", "a") as file:
        file.write("""def new_compute_dHdp1(m1, m2, EMgamma, tortoise, x, y, z, p1, p2, p3, S1x, S1y, S1z, S2x, S2y, S2z):
""")
        for i in range(len(lr) - 1):
            file.write("    " + lr[i].lhs + " = " + str(lr[i].rhs).replace("Rational(", "np.true_divide(").replace("sqrt(", "np.sqrt(").replace("log(", "np.log(").replace("sign(", "np.sign(").replace("Abs(", "np.abs(").replace("pi", "np.pi") + "\n")
        for i in range(len(lhss_deriv_p1)):
            file.write("    " + str(lhss_deriv_p1[i]).replace("prm", "prm_p1") + " = " + replace_numpy_funcs(rhss_deriv_p1[i]).replace("prm", "prm_p1").replace("sp.sqrt(","np.sqrt(").replace("sp.log(","np.log(").replace("sp.sign(","np.sign(").replace("sp.Abs(", "np.abs(") + "\n")
        file.write("    return np.array([Hrealprm_p1])\n")        
        
    with open("SEOBNR_Playground_Pycodes/new_dH.py", "a") as file:
        file.write("""def new_compute_dHdp2(m1, m2, EMgamma, tortoise, x, y, z, p1, p2, p3, S1x, S1y, S1z, S2x, S2y, S2z):
""")
        for i in range(len(lr) - 1):
            file.write("    " + lr[i].lhs + " = " + str(lr[i].rhs).replace("Rational(", "np.true_divide(").replace("sqrt(", "np.sqrt(").replace("log(", "np.log(").replace("sign(", "np.sign(").replace("Abs(", "np.abs(").replace("pi", "np.pi") + "\n")
        for i in range(len(lhss_deriv_p2)):
            file.write("    " + str(lhss_deriv_p2[i]).replace("prm", "prm_p2") + " = " + replace_numpy_funcs(rhss_deriv_p2[i]).replace("prm", "prm_p2").replace("sp.sqrt(","np.sqrt(").replace("sp.log(","np.log(").replace("sp.sign(","np.sign(").replace("sp.Abs(", "np.abs(") + "\n")
        file.write("    return np.array([Hrealprm_p2])\n")

    with open("SEOBNR_Playground_Pycodes/new_dH.py", "a") as file:
        file.write("""def new_compute_dHdp3(m1, m2, EMgamma, tortoise, x, y, z, p1, p2, p3, S1x, S1y, S1z, S2x, S2y, S2z):
""")
        for i in range(len(lr) - 1):
            file.write("    " + lr[i].lhs + " = " + str(lr[i].rhs).replace("Rational(", "np.true_divide(").replace("sqrt(", "np.sqrt(").replace("log(", "np.log(").replace("sign(", "np.sign(").replace("Abs(", "np.abs(").replace("pi", "np.pi") + "\n")
        for i in range(len(lhss_deriv_p3)):
            file.write("    " + str(lhss_deriv_p3[i]).replace("prm", "prm_p3") + " = " + replace_numpy_funcs(rhss_deriv_p3[i]).replace("prm", "prm_p3").replace("sp.sqrt(","np.sqrt(").replace("sp.log(","np.log(").replace("sp.sign(","np.sign(").replace("sp.Abs(", "np.abs(") + "\n")
        file.write("    return np.array([Hrealprm_p3])\n")
        
    with open("SEOBNR_Playground_Pycodes/new_dH.py", "a") as file:
        file.write("""def new_compute_dHdS1x(m1, m2, EMgamma, tortoise, x, y, z, p1, p2, p3, S1x, S1y, S1z, S2x, S2y, S2z):
""")
        for i in range(len(lr) - 1):
            file.write("    " + lr[i].lhs + " = " + str(lr[i].rhs).replace("Rational(", "np.true_divide(").replace("sqrt(", "np.sqrt(").replace("log(", "np.log(").replace("sign(", "np.sign(").replace("Abs(", "np.abs(").replace("pi", "np.pi") + "\n")
        for i in range(len(lhss_deriv_S1x)):
            file.write("    " + str(lhss_deriv_S1x[i]).replace("prm", "prm_S1x") + " = " + replace_numpy_funcs(rhss_deriv_S1x[i]).replace("prm", "prm_S1x").replace("sp.sqrt(","np.sqrt(").replace("sp.log(","np.log(").replace("sp.sign(","np.sign(").replace("sp.Abs(", "np.abs(") + "\n")
        file.write("    return np.array([Hrealprm_S1x])\n")
    
    with open("SEOBNR_Playground_Pycodes/new_dH.py", "a") as file:
        file.write("""def new_compute_dHdS1y(m1, m2, EMgamma, tortoise, x, y, z, p1, p2, p3, S1x, S1y, S1z, S2x, S2y, S2z):
""")
        for i in range(len(lr) - 1):
            file.write("    " + lr[i].lhs + " = " + str(lr[i].rhs).replace("Rational(", "np.true_divide(").replace("sqrt(", "np.sqrt(").replace("log(", "np.log(").replace("sign(", "np.sign(").replace("Abs(", "np.abs(").replace("pi", "np.pi") + "\n")
        for i in range(len(lhss_deriv_S1y)):
            file.write("    " + str(lhss_deriv_S1y[i]).replace("prm", "prm_S1y") + " = " + replace_numpy_funcs(rhss_deriv_S1y[i]).replace("prm", "prm_S1y").replace("sp.sqrt(","np.sqrt(").replace("sp.log(","np.log(").replace("sp.sign(","np.sign(").replace("sp.Abs(", "np.abs(") + "\n")
        file.write("    return np.array([Hrealprm_S1y])\n")
    
    with open("SEOBNR_Playground_Pycodes/new_dH.py", "a") as file:
        file.write("""def new_compute_dHdS1z(m1, m2, EMgamma, tortoise, x, y, z, p1, p2, p3, S1x, S1y, S1z, S2x, S2y, S2z):
""")
        for i in range(len(lr) - 1):
            file.write("    " + lr[i].lhs + " = " + str(lr[i].rhs).replace("Rational(", "np.true_divide(").replace("sqrt(", "np.sqrt(").replace("log(", "np.log(").replace("sign(", "np.sign(").replace("Abs(", "np.abs(").replace("pi", "np.pi") + "\n")
        for i in range(len(lhss_deriv_S1z)):
            file.write("    " + str(lhss_deriv_S1z[i]).replace("prm", "prm_S1z") + " = " + replace_numpy_funcs(rhss_deriv_S1z[i]).replace("prm", "prm_S1z").replace("sp.sqrt(","np.sqrt(").replace("sp.log(","np.log(").replace("sp.sign(","np.sign(").replace("sp.Abs(", "np.abs(") + "\n")
        file.write("    return np.array([Hrealprm_S1z])\n")
    
    with open("SEOBNR_Playground_Pycodes/new_dH.py", "a") as file:
        file.write("""def new_compute_dHdS2x(m1, m2, EMgamma, tortoise, x, y, z, p1, p2, p3, S1x, S1y, S1z, S2x, S2y, S2z):
""")
        for i in range(len(lr) - 1):
            file.write("    " + lr[i].lhs + " = " + str(lr[i].rhs).replace("Rational(", "np.true_divide(").replace("sqrt(", "np.sqrt(").replace("log(", "np.log(").replace("sign(", "np.sign(").replace("Abs(", "np.abs(").replace("pi", "np.pi") + "\n")
        for i in range(len(lhss_deriv_S2x)):
            file.write("    " + str(lhss_deriv_S2x[i]).replace("prm", "prm_S2x") + " = " + replace_numpy_funcs(rhss_deriv_S2x[i]).replace("prm", "prm_S2x").replace("sp.sqrt(","np.sqrt(").replace("sp.log(","np.log(").replace("sp.sign(","np.sign(").replace("sp.Abs(", "np.abs(") + "\n")
        file.write("    return np.array([Hrealprm_S2x])\n")
    
    with open("SEOBNR_Playground_Pycodes/new_dH.py", "a") as file:
        file.write("""def new_compute_dHdS2y(m1, m2, EMgamma, tortoise, x, y, z, p1, p2, p3, S1x, S1y, S1z, S2x, S2y, S2z):
""")
        for i in range(len(lr) - 1):
            file.write("    " + lr[i].lhs + " = " + str(lr[i].rhs).replace("Rational(", "np.true_divide(").replace("sqrt(", "np.sqrt(").replace("log(", "np.log(").replace("sign(", "np.sign(").replace("Abs(", "np.abs(").replace("pi", "np.pi") + "\n")
        for i in range(len(lhss_deriv_S2y)):
            file.write("    " + str(lhss_deriv_S2y[i]).replace("prm", "prm_S2y") + " = " + replace_numpy_funcs(rhss_deriv_S2y[i]).replace("prm", "prm_S2y").replace("sp.sqrt(","np.sqrt(").replace("sp.log(","np.log(").replace("sp.sign(","np.sign(").replace("sp.Abs(", "np.abs(") + "\n")
        file.write("    return np.array([Hrealprm_S2y])\n")
    
    with open("SEOBNR_Playground_Pycodes/new_dH.py", "a") as file:
        file.write("""def new_compute_dHdS2z(m1, m2, EMgamma, tortoise, x, y, z, p1, p2, p3, S1x, S1y, S1z, S2x, S2y, S2z):
""")
        for i in range(len(lr) - 1):
            file.write("    " + lr[i].lhs + " = " + str(lr[i].rhs).replace("Rational(", "np.true_divide(").replace("sqrt(", "np.sqrt(").replace("log(", "np.log(").replace("sign(", "np.sign(").replace("Abs(", "np.abs(").replace("pi", "np.pi") + "\n")
        for i in range(len(lhss_deriv_S2z)):
            file.write("    " + str(lhss_deriv_S2z[i]).replace("prm", "prm_S2z") + " = " + replace_numpy_funcs(rhss_deriv_S2z[i]).replace("prm", "prm_S2z").replace("sp.sqrt(","np.sqrt(").replace("sp.log(","np.log(").replace("sp.sign(","np.sign(").replace("sp.Abs(", "np.abs(") + "\n")
        file.write("    return np.array([Hrealprm_S2z])\n")

### create a functions that calls and returns all derivatives for integrations

    with open("SEOBNR_Playground_Pycodes/new_dH.py", "a") as file:
        file.write("""def compute_Hreal_and_csi_derivatives(m1, m2, EMgamma, tortoise, x, y, z, p1, p2, p3, S1x, S1y, S1z, S2x, S2y, S2z):
    dHdx = new_compute_dHdx(m1, m2, EMgamma, tortoise, x, y, z, p1, p2, p3, S1x, S1y, S1z, S2x, S2y, S2z)
    dHdy = new_compute_dHdy(m1, m2, EMgamma, tortoise, x, y, z, p1, p2, p3, S1x, S1y, S1z, S2x, S2y, S2z)
    dHdz = new_compute_dHdz(m1, m2, EMgamma, tortoise, x, y, z, p1, p2, p3, S1x, S1y, S1z, S2x, S2y, S2z)
    dHdpx = new_compute_dHdp1(m1, m2, EMgamma, tortoise, x, y, z, p1, p2, p3, S1x, S1y, S1z, S2x, S2y, S2z)
    dHdpy = new_compute_dHdp2(m1, m2, EMgamma, tortoise, x, y, z, p1, p2, p3, S1x, S1y, S1z, S2x, S2y, S2z)
    dHdpz = new_compute_dHdp3(m1, m2, EMgamma, tortoise, x, y, z, p1, p2, p3, S1x, S1y, S1z, S2x, S2y, S2z)
    dHdS1x = new_compute_dHdS1x(m1, m2, EMgamma, tortoise, x, y, z, p1, p2, p3, S1x, S1y, S1z, S2x, S2y, S2z)
    dHdS1y = new_compute_dHdS1y(m1, m2, EMgamma, tortoise, x, y, z, p1, p2, p3, S1x, S1y, S1z, S2x, S2y, S2z)
    dHdS1z = new_compute_dHdS1z(m1, m2, EMgamma, tortoise, x, y, z, p1, p2, p3, S1x, S1y, S1z, S2x, S2y, S2z)
    dHdS2x = new_compute_dHdS2x(m1, m2, EMgamma, tortoise, x, y, z, p1, p2, p3, S1x, S1y, S1z, S2x, S2y, S2z)
    dHdS2y = new_compute_dHdS2y(m1, m2, EMgamma, tortoise, x, y, z, p1, p2, p3, S1x, S1y, S1z, S2x, S2y, S2z)
    dHdS2z = new_compute_dHdS2z(m1, m2, EMgamma, tortoise, x, y, z, p1, p2, p3, S1x, S1y, S1z, S2x, S2y, S2z)
    return np.array([dHdx[0],dHdy[0],dHdz[0],dHdpx[0],dHdpy[0],dHdpz[0],dHdS1x[0],dHdS1y[0],dHdS1z[0],dHdS2x[0],dHdS2y[0],dHdS2z[0],dHdx[1],dHdy[1],dHdz[1]])
""")
    
    # TylerK: now create a text file listing only the terms so we can take a second derivative!

    with open("SEOBNR_Playground_Pycodes/dHdx.txt", "w") as file:
        for i in range(len(lr) - 1):
            file.write(lr[i].lhs + " = " + str(lr[i].rhs) + "\n")
        for i in range(len(lhss_deriv_x)):
            file.write(str(lhss_deriv_x[i]).replace("prm", "prm_x") + " = " + replace_numpy_funcs(rhss_deriv_x[i]).replace("prm", "prm_x") + "\n")
    
    with open("SEOBNR_Playground_Pycodes/dHdp2.txt", "w") as file:
        for i in range(len(lr) - 1):
            file.write(lr[i].lhs + " = " + str(lr[i].rhs) + "\n")
        for i in range(len(lhss_deriv_p2)):
            file.write(str(lhss_deriv_p2[i]).replace("prm", "prm_p2") + " = " + replace_numpy_funcs(rhss_deriv_p2[i]).replace("prm", "prm_p2") + "\n")

    with open("SEOBNR_Playground_Pycodes/dHdp3.txt", "w") as file:
        for i in range(len(lr) - 1):
            file.write(lr[i].lhs + " = " + str(lr[i].rhs) + "\n")
        for i in range(len(lhss_deriv_p3)):
            file.write(str(lhss_deriv_p3[i]).replace("prm", "prm_p3") + " = " + replace_numpy_funcs(rhss_deriv_p3[i]).replace("prm", "prm_p3") + "\n")
