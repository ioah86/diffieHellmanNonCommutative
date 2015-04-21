"""
This is a Sage script that allows the user to generate
automatically the header and the source C-file for a certain finite
field. This finite field has to fulfill the property that it has at
least size p^k, where p is a prime number and k >1 is an integer.
The generated C-files can directly be plugged in into the
Ore Diffie-Hellman implementation.

:Author: Albert Heinle <albert.heinle@googlemail.com>
"""

f = open("gf_coefficients_h_template.txt","r")
gf_coefficients_h_template = f.read()
f.close()

f = open("gf_coefficients_c_template.txt","r")
gf_coefficients_c_template = f.read()
f.close()

def is_proper_gf(f):
    """
    ring -> bool
    This function takes a field as an input (f), and returns true if
    all of the following conditions are fulfilled:
    - F is a finite field
    - F has p^k elements, where p is the characteristic and k is an
      integer greater than 1
    """
    if not f.is_field():
        return False
    if not f.is_finite():
        return False
    if f.is_prime_field():
        return False
    return True

def make_gf_coefficients_h(f):
    """
    ring -> string
    This functions takes a field f as an input, and begins with making
    certain sanity checks on it. If it fulfills all requirements to be
    used in our context (i.e. is_proper_gf(f) returns True), this
    function will create the corresponding header-file for it to serve
    as c-code.
    """
    if not is_proper_gf(f):
        return ""
    modulus = f.characteristic()
    degree_of_extension = f.degree()
    number_of_elements_in_f = f.order()
    homs = Hom(f,f)
    homs_header_lines = ""
    nhom = 1
    for h in homs:
        if h.is_identity():
            continue
        hTemp = h
        hOrder = 1
        while not hTemp.is_identity():
            hTemp = hTemp * h
            hOrder = hOrder +1
        homs_header_lines += """struct GFModulus Hom%i(struct GFModulus);
#define ORDERHOM%i (%i)
"""%(nhom,nhom,hOrder)
        nhom = nhom +1
    result = gf_coefficients_h_template % (modulus,
                                           degree_of_extension,
                                           number_of_elements_in_f,
                                           homs_header_lines)
    return result

def make_gf_coefficients_c(f):
    """
    ring -> string
    This functions takes a field f as an input, and begins with making
    certain sanity checks on it. If it fulfills all requirements to be
    used in our context (i.e. is_proper_gf(f) returns True), this
    function will create the corresponding c-file for its elements and
    arithmetics.
    """
    if not is_proper_gf(f):
        return ""
    #getting minimal polynomial
    minPolyInString = str(F.polynomial())
    #getting the body of the function for multiplying two polynomials
    a = [var('a_%d'%i) for i in range(f.degree())]
    b = [var('b_%d'%i) for i in range(f.degree())]
    multGFbody = "";
    for i in range(len(a)):
        multGFbody+=("  int %s = inp1.coeffs[%i];\n"%(str(a[i]),i))
    for i in range(len(b)):
        multGFbody+=("  int %s = inp2.coeffs[%i];\n"%(str(b[i]),i))
    factor1 = sum([a[i]*f.gen()^i for i in range(f.degree())])
    factor2 = sum([b[i]*f.gen()^i for i in range(f.degree())])
    multRes = (expand(factor1*factor2)).coefficients(f.gen())
    for i in multRes:
        multGFbody +="  result.coeffs[%i] = (%s) \
%%MODULUS;\n"%(i[1],str(i[0]))
    #Dealing with the homs
    homs = Hom(f,f)
    homs_code_blocks = ""
    nhom = 1
    for h in homs:
        if h.is_identity():
            continue
        hTemp = h
        homs_code_blocks += """struct GFModulus Hom%i(struct GFModulus inp)
{
  //Short description of Hom%i:
  //%s |--> %s
  struct GFModulus result;
"""%(nhom,nhom,str(f.gen()), str(h(f.gen())))
        for i in range(len(a)):
            homs_code_blocks+=("  int %s = \
inp.coeffs[%i];\n"%(str(a[i]),i))
        factorTemp = (sum([a[i]*h(f.gen())^i for i in
                               range(f.degree())])).coefficients(f.gen())
        for i in factorTemp:
            homs_code_blocks+="  result.coeffs[%i] = (%s) \
%%MODULUS;\n"%(i[1],str(i[0]))
        homs_code_blocks+="""  return result;
}

"""
        nhom = nhom +1
    result = gf_coefficients_c_template%(minPolyInString, multGFbody,
                                         homs_code_blocks)
    return result
