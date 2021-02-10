# NRPy+ code to generate second derivatives of the SEOBNRv3_Sph Hamiltonian from a list of numerical expressions computing
# said Hamiltonian. Originally written by Zach Etienne; edited and commented by Tyler Knowles; Adapted to spherical by Sid Mahesh

from outputC import outputC,lhrh,superfast_uniq  # NRPy+: Core C code output module
import sympy as sp      # SymPy: The Python computer algebra package upon which NRPy+ depends
import sys              # Python module for multiplatform OS-related functions


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
# example, pass 'rprm=1' to replace each instance of 'rprm' with 1 and 'qprm' with 0 for each q in (theta,phi,pr,ptheta,pphi,S1r,
# S1theta,S1phi,S2r,S2theta,S2phi).  This produces expressions which compute the partial derivative of the Hamiltonian with respect
# to r.
def deriv_onevar(lhss_deriv, rhss_deriv, rprm=0, thetaprm=0, phiprm=0, prprm=0, pthetaprm=0, pphiprm=0, S1rprm=0, S1thetaprm=0,
                 S1phiprm=0, S2rprm=0, S2thetaprm=0, S2phiprm=0):
    if rprm + thetaprm + phiprm + prprm + pthetaprm + pphiprm + S1rprm + S1thetaprm + S1phiprm + S2rprm + S2thetaprm + S2phiprm != 1:
        print("deriv_onevar() cannot take more than one derivative at a time!")
        sys.exit()

    # Create 'new' arrays to store and manipulate derivative terms.
    lhss_deriv_new = []
    rhss_deriv_new = []
    # Append derivative terms to 'new' arrays
    for i in range(len(rhss_deriv)):
        lhss_deriv_new.append(lhss_deriv[i])
        rhss_deriv_new.append(rhss_deriv[i])
    # Replace each instance of 'qprm', q in (r,theta,phi,pr,ptheta,pphi,S1r,S1theta,S1phi,S2r,S2theta,S2phi), with either 0 or 1.
    for i in range(len(rhss_deriv_new)):
        for var in rhss_deriv_new[i].free_symbols:
            if str(var) == "rprm":
                rhss_deriv_new[i] = rhss_deriv_new[i].subs(var, rprm)
            elif str(var) == "thetaprm":
                rhss_deriv_new[i] = rhss_deriv_new[i].subs(var, thetaprm)
            elif str(var) == "phiprm":
                rhss_deriv_new[i] = rhss_deriv_new[i].subs(var, phiprm)
            elif str(var) == "prprm":
                rhss_deriv_new[i] = rhss_deriv_new[i].subs(var, prprm)
            elif str(var) == "pthetaprm":
                rhss_deriv_new[i] = rhss_deriv_new[i].subs(var, pthetaprm)
            elif str(var) == "pphiprm":
                rhss_deriv_new[i] = rhss_deriv_new[i].subs(var, pphiprm)
            elif str(var) == "S1rprm":
                rhss_deriv_new[i] = rhss_deriv_new[i].subs(var, S1rprm)
            elif str(var) == "S1thetaprm":
                rhss_deriv_new[i] = rhss_deriv_new[i].subs(var, S1thetaprm)
            elif str(var) == "S1phiprm":
                rhss_deriv_new[i] = rhss_deriv_new[i].subs(var, S1phiprm)
            elif str(var) == "S2rprm":
                rhss_deriv_new[i] = rhss_deriv_new[i].subs(var, S2rprm)
            elif str(var) == "S2thetaprm":
                rhss_deriv_new[i] = rhss_deriv_new[i].subs(var, S2thetaprm)
            elif str(var) == "S2phiprm":
                rhss_deriv_new[i] = rhss_deriv_new[i].subs(var, S2phiprm)
    # Simplify the derivative expressions with simplify_deriv().
    lhss_deriv_simp, rhss_deriv_simp = simplify_deriv(lhss_deriv_new, rhss_deriv_new)
    # Return simplified derivative expression.
    return lhss_deriv_simp, rhss_deriv_simp
    #return lhss_deriv_new, rhss_deriv_new


# replace_numpy_funcs() replaces specific SymPy function names with the corresponding NumPy function names.
def replace_numpy_funcs(expression):
    return str(expression).replace("sqrt(", "sp.sqrt(").replace("Abs(", "sp.Abs(").replace("log(",
                                                                "sp.log(").replace("sign(", "sp.sign(").replace("sin(","sp.sin(")


# output_H_sec_derivs() is the main wrapper function for computing the SEONBRv3_Sph Hamiltonian H and the second
# partial derivatives of H for r, r and r , pphi for initial conditions

def output_H_sec_derivs():
    # Open and read the file of numerical expressions (written in SymPy syntax) computing the SEOBNRv3_Sph Hamiltonian derivative.
    
    f = open("Playground_Pycodes/dHdr.txt", 'r')
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
            lr.append(lhrh(lhs=splitHamterms[0].replace(" ", "").replace("Lambda", "Lamb").replace("prm", ""),
                           rhs=splitHamterms[1].replace(" ", "").replace("sp.", "").replace("Lambda", "Lamb").replace("prm", "")))
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
        lhss.append(sp.sympify(lr[i].lhs))
        rhss.append(sp.sympify(lr[i].rhs))
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
        lhss_deriv.append(sp.sympify(str(lhss[i]) + "prm"))
        newrhs = sp.sympify(
            str(sp.diff(rhss[i], xx)).replace("(xx)", "").replace(", xx", "prm").replace("Derivative", ""))
        rhss_deriv.append(newrhs)
    # Simplify derivative expressions with simplify_deriv()
    lhss_deriv_simp, rhss_deriv_simp = simplify_deriv(lhss_deriv, rhss_deriv)
    lhss_deriv = lhss_deriv_simp
    rhss_deriv = rhss_deriv_simp
    # Generate partial derivatives with respect to r and pphi only
    lhss_deriv_r, rhss_deriv_r = deriv_onevar(lhss_deriv, rhss_deriv, rprm=1, thetaprm=0, phiprm=0, prprm=0, pthetaprm=0, pphiprm=0,
                                              S1rprm=0, S1thetaprm=0, S1phiprm=0, S2rprm=0, S2thetaprm=0, S2phiprm=0)
    lhss_deriv_pphi, rhss_deriv_pphi = deriv_onevar(lhss_deriv, rhss_deriv, rprm=0, thetaprm=0, phiprm=0, prprm=0, pthetaprm=0, pphiprm=1,
                                              S1rprm=0, S1thetaprm=0, S1phiprm=0, S2rprm=0, S2thetaprm=0, S2phiprm=0)
    
    # Prepare to output derivative expressions in C syntax
#    outstring = "/* SEOBNR Hamiltonian expression: */\n"
#    outstringsp = ""
#    outsplhs = []
#    outsprhs = []
#    for i in range(len(lr)):
#        outstring += outputC(sp.sympify(lr[i].rhs), lr[i].lhs, "returnstring",
#                             "outCverbose=False,includebraces=False,CSE_enable=False")
#        outstringsp += lr[i].lhs + " = " + lr[i].rhs + "\n"
#        outsplhs.append(sp.sympify(lr[i].lhs))
#        outsprhs.append(sp.sympify(lr[i].rhs))

#    outstring += "\n\n\n/* SEOBNR \partial_x H expression: */\n"
#    for i in range(len(lhss_deriv_x)):
#        outstring += outputC(rhss_deriv_x[i], str(lhss_deriv_x[i]), "returnstring",
#                             "outCverbose=False,includebraces=False,CSE_enable=False")
#        outstringsp += str(lhss_deriv_x[i]) + " = " + str(rhss_deriv_x[i]) + "\n"
#        outsplhs.append(lhss_deriv_x[i])
#        outsprhs.append(rhss_deriv_x[i])

    with open("Playground_Pycodes/d2Hdr2.py", "w") as file:
        file.write("""from __future__ import division
import numpy as np
def compute_d2Hdr2(m1, m2, eta, r, theta, phi, pr, ptheta, pphi, S1r, S1theta, S1phi, S2r, S2theta, S2phi, KK, k0, k1, dSO, dSS, tortoise, EMgamma):
""")
        for i in range(len(lr)):
            file.write("    " + lr[i].lhs + " = " + str(lr[i].rhs).replace("Rational(", "np.true_divide(").replace("sqrt(", "np.sqrt(").replace("log(", "np.log(").replace("sign(", "np.sign(").replace("Abs(", "np.abs(").replace("pi", "np.pi").replace("sin(","np.sin(") + "\n")
        for i in range(len(lhss_deriv_r)):
            file.write("    " + str(lhss_deriv_r[i]).replace("prm", "prm_r") + " = " + replace_numpy_funcs(rhss_deriv_r[i]).replace("prm", "prm_r").replace("sp.sqrt(","np.sqrt(").replace("sp.log(","np.log(").replace("sp.sign(","np.sign(").replace("sp.Abs(", "np.abs(").replace("sp.sin(","np.sin(") + "\n")
        #file.write("    Hreal=np.sqrt(1+2*eta*(Heff-1))"+"\n")    
        file.write("    return np.array([Hreal,Hreal_rprm_r])")    

    with open("Playground_Pycodes/d2Hdrdpphi.py", "w") as file:
        file.write("""from __future__ import division
import numpy as np
def compute_d2Hdrdpphi(m1, m2, eta, r, theta, phi, pr, ptheta, pphi, S1r, S1theta, S1phi, S2r, S2theta, S2phi, KK, k0, k1, dSO, dSS, tortoise, EMgamma):
""")
        for i in range(len(lr) - 1):
            file.write("    " + lr[i].lhs + " = " + str(lr[i].rhs).replace("Rational(", "np.true_divide(").replace("sqrt(", "np.sqrt(").replace("log(", "np.log(").replace("sign(", "np.sign(").replace("Abs(", "np.abs(").replace("pi", "np.pi").replace("sin(","np.sin(") + "\n")
        for i in range(len(lhss_deriv_pphi)):
            file.write("    " + str(lhss_deriv_pphi[i]).replace("prm", "prm_pphi") + " = " + replace_numpy_funcs(rhss_deriv_pphi[i]).replace("prm", "prm_pphi").replace("sp.sqrt(","np.sqrt(").replace("sp.log(","np.log(").replace("sp.sign(","np.sign(").replace("sp.Abs(", "np.abs(").replace("sp.sin(","np.sin(") + "\n")
        file.write("    return np.array([Hreal_rprm_pphi])")
    
    f = open("Playground_Pycodes/dHdpphi.txt", 'r')
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
            lr.append(lhrh(lhs=splitHamterms[0].replace(" ", "").replace("Lambda", "Lamb").replace("prm", ""),
                           rhs=splitHamterms[1].replace(" ", "").replace("sp.", "").replace("Lambda", "Lamb").replace("prm", "")))
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
        lhss.append(sp.sympify(lr[i].lhs))
        rhss.append(sp.sympify(lr[i].rhs))
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
        lhss_deriv.append(sp.sympify(str(lhss[i]) + "prm"))
        newrhs = sp.sympify(
            str(sp.diff(rhss[i], xx)).replace("(xx)", "").replace(", xx", "prm").replace("Derivative", ""))
        rhss_deriv.append(newrhs)
    # Simplify derivative expressions with simplify_deriv()
    lhss_deriv_simp, rhss_deriv_simp = simplify_deriv(lhss_deriv, rhss_deriv)
    lhss_deriv = lhss_deriv_simp
    rhss_deriv = rhss_deriv_simp
    # Generate partial derivatives with respect to each of the twelve input variables
    lhss_deriv_r, rhss_deriv_r = deriv_onevar(lhss_deriv, rhss_deriv, rprm=1, thetaprm=0, phiprm=0, prprm=0, pthetaprm=0, pphiprm=0,
                                              S1rprm=0, S1thetaprm=0, S1phiprm=0, S2rprm=0, S2thetaprm=0, S2phiprm=0)
    # Prepare to output derivative expressions in C syntax
#    outstring = "/* SEOBNR Hamiltonian expression: */\n"
#    outstringsp = ""
#    outsplhs = []
#    outsprhs = []
#    for i in range(len(lr)):
#        outstring += outputC(sp.sympify(lr[i].rhs), lr[i].lhs, "returnstring",
#                             "outCverbose=False,includebraces=False,CSE_enable=False")
#        outstringsp += lr[i].lhs + " = " + lr[i].rhs + "\n"
#        outsplhs.append(sp.sympify(lr[i].lhs))
#        outsprhs.append(sp.sympify(lr[i].rhs))
#
#    outstring += "\n\n\n/* SEOBNR \partial_x H expression: */\n"
#    for i in range(len(lhss_deriv_x)):
#        outstring += outputC(rhss_deriv_x[i], str(lhss_deriv_x[i]), "returnstring",
#                             "outCverbose=False,includebraces=False,CSE_enable=False")
#        outstringsp += str(lhss_deriv_x[i]) + " = " + str(rhss_deriv_x[i]) + "\n"
#        outsplhs.append(lhss_deriv_x[i])
#        outsprhs.append(rhss_deriv_x[i])
#
    with open("Playground_Pycodes/d2Hdpphidr.py", "w") as file:
        file.write("""from __future__ import division
import numpy as np
def compute_d2Hdpphidr(m1, m2, eta, r, theta, phi, pr, ptheta, pphi, S1r, S1theta, S1phi, S2r, S2theta, S2phi, KK, k0, k1, dSO, dSS, tortoise, EMgamma):
""")
        for i in range(len(lr) - 1):
            file.write("    " + lr[i].lhs + " = " + str(lr[i].rhs).replace("Rational(", "np.true_divide(").replace("sqrt(", "np.sqrt(").replace("log(", "np.log(").replace("sign(", "np.sign(").replace("Abs(", "np.abs(").replace("pi", "np.pi").replace("sin(","np.sin(") + "\n")
        for i in range(len(lhss_deriv_r)):
            file.write("    " + str(lhss_deriv_r[i]).replace("prm", "prm_r") + " = " + replace_numpy_funcs(rhss_deriv_r[i]).replace("prm", "prm_r").replace("sp.sqrt(","np.sqrt(").replace("sp.log(","np.log(").replace("sp.sign(","np.sign(").replace("sp.Abs(", "np.abs(").replace("sp.sin(","np.sin(") + "\n")
        file.write("    return np.array([Hreal_pphiprm_r])")    
