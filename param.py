# Given a parameter definition, generates necessary C++, Python and XML bindings
# coding=utf-8

def main():
    # read the input definition
    # name, description, datatype, default
    PARAMS = {
        'common': [
            ('verbose', '', 'bool', False),
#            ('verbose', '', 'bool', True),
            ('output_path', 'directory where model outputs should be placed (must end in /)', 'string', 'output-Uranus/'),

#            ('paraview_flag','flag to control if create paraview panorama', 'bool', False),
            ('paraview_flag','flag to control if create paraview panorama', 'bool', True),
       ],






        'uranus': [

#            ('nm', 'the maximum number of iterations', 'int', 4),
            ('nm', 'the maximum number of iterations', 'int', 224),
#            ('nm', 'the maximum number of iterations', 'int', 512),
#            ('checkpoint', "control when to write output files", 'int', 32),
#            ('checkpoint', "control when to write output files", 'int', 16),
            ('checkpoint', "control when to write output files", 'int', 8),
#            ('checkpoint', "control when to write output files", 'int', 2),

            ('panorama_print', "control when to write panorama files", 'int', 32),
#            ('panorama_print', "control when to write panorama files", 'int', 256),


            ('Coriolis', 'Coriolis force', 'double', 1),
            ('centrifugal', 'centrifugal force', 'double', 1),
            ('buoyancy', 'buoyancy force', 'double', 1),
#            ('buoyancy', 'buoyancy force', 'double', 0),
            ('chemical_reaction', 'chemical reactions included', 'double', 1),
            ('epsres', 'accuracy of relative and absolute errors', 'double', 0.00001),

            ('L_atm', 'extension of the troposhere in km, 360km/40 steps = 9km', 'double', 360.0),

#            ('tropopause_pole', 'extension of the troposphere at the poles in km', 'double', 115.0),
#            ('tropopause_equator', 'extension of the troposphere at the equator in km', 'double', 125.0),

            ('tropopause_pole', 'extension of the troposphere at the poles in km', 'double', 125.0),
            ('tropopause_equator', 'extension of the troposphere at the equator in km', 'double', 115.0),

            ('re', 'Reynolds number: ratio viscous to inertia forces, Re = u * L/nue', 'double', 1000.0),
            ('ec', 'Eckert number: ratio kinetic energy to enthalpy, Ec = u²/cp T', 'double', 0.00044),

            ('ep', 'ratio of the gas constants of dry air to water vapour [/]', 'double', 0.623),
            ('hp', 'water vapour pressure at T = 0°C: E = 6.1 hPa', 'double', 6.1078),
            ('lv', 'specific latent evaporation heat(condensation heat) in J/kg', 'double', 2.52e6),
            ('ls', 'specific latent vaporisation heat(sublimation heat) in J/kg', 'double', 2.83e6),
            ('cp_l', 'specific heat capacity of dry air at constant pressure and 20°C in J/(kg K)', 'double', 1005.0),

            ('pr', 'Prandtl number of h2oe for laminar flows', 'double', 0.69),
            ('g', 'gravitational acceleration of Uranus in m/s²', 'double', 8.69),
            ('omega', 'rotation number of Uranus in 1/s', 'double', 1.01e-4),

            ('u_0', 'average value of velocity in m/s', 'double', 100.0),
            ('r_0_water', 'reference density of fresh water in kg/m3', 'double', 997.0),

            ('ua', 'initial velocity component in r-direction', 'double', 0.0),
            ('va', 'initial velocity component in theta-direction', 'double', 0.0),
            ('wa', 'initial velocity component in phi-direction', 'double', 0.0),
            ('pa', 'initial value for the pressure field', 'double', 0.0),
            ('ca', 'value 0.04 stands for the maximum value of 40 g/kg water vapour', 'double', 0.0),
            ('ta', 'initial value for the temperature field, 1.0 compares to 0° C compares to 273.15 K', 'double', 1.0),

            ('t_ref', 'temperature in K compare to -196.75°C', 'double', 76.4),
            ('t_equator', 'temperature at the equator -212.5°C compares to 60.65K', 'double', 60.65),
            ('t_pole', 'temperature at the poles - assumption -224°C compares to 49.15K', 'double', 49.15),

            ('p_ref', 'pressure in bar', 'double', 1.0),

            ('R_ref', 'average gas constant in J/(g*K) after Sanchez et. al.', 'double', 3.61491),

            ('ch4_tropopause', 'minimum ch4 at tropopause kg/kg', 'double', 7e-5),
            ('h2o_tropopause', 'minimum h2o vapour at tropopause kg/kg', 'double', 6.9e-10),
            ('h2s_tropopause', 'minimum h2s at tropopause kg/kg', 'double', 1.5e-9),
            ('nh3_tropopause', 'minimum rate nh3 at tropopause kg/kg', 'double', 7.4e-13),
            ('nh4sh_tropopause', 'salt nh4sh at tropopause kg/kg', 'double', 0.0),

       ],
 
    }

    XML_READ_FUNCS = {
        "string": "FillStringWithElement",
        "double": "FillDoubleWithElement",
        "int": "FillIntWithElement",
        "bool": "FillBoolWithElement"
    }

    def write_cpp_defaults(filename, classname, sections):
        with open(filename, 'w') as f:
            f.write("// header files\n")
            f.write("// THIS FILE IS AUTOMATICALLY GENERATED BY param.py\n")
            f.write("// ANY CHANGES WILL BE OVERWRITTEN AT COMPILE TIME\n")
            f.write("\n")
            f.write("void %s::SetDefaultConfig() {\n" % classname)
            for section in sections:
                f.write('\n  // %s section\n' % section)
                for slug, desc, ctype, default in PARAMS[section]:
                    rhs = default
                    if ctype == 'string':
                        rhs = '"%s"' % default
                    elif ctype == 'bool':
                        if default:
                            rhs = 'true'
                        else:
                            rhs = 'false'
                    f.write('  %s = %s;\n' % (slug, rhs))
            f.write("}")

    def write_cpp_load_config(filename, classname, sections):
        with open(filename, 'w') as f:
            f.write("// config files\n")
            f.write("// THIS FILE IS AUTOMATICALLY GENERATED BY param.py\n")
            f.write("// ANY CHANGES WILL BE OVERWRITTEN AT COMPILE TIME\n")
            f.write("\n")
            for section in sections:
                f.write('\n  // %s section\n' % section)
                element_var_name = 'elem_%s' % section
                f.write('\n  if (%s) {\n' % (element_var_name))
                for slug, desc, ctype, default in PARAMS [section]:
                    func_name = XML_READ_FUNCS [ctype]
                    f.write('    Config::%s(%s, "%s", %s);\n' % (func_name, element_var_name, slug, slug))
                f.write("  }\n")

    def write_cpp_ueaders(filename, sections, is_extern=False):
        with open(filename, 'w') as f:
            f.write("// header files\n")
            f.write("// THIS FILE IS AUTOMATICALLY GENERATED BY param.py\n")
            f.write("// ANY CHANGES WILL BE OVERWRITTEN AT COMPILE TIME\n")
            f.write("\n")
            if is_extern:
                f.write("#include<string>\n\n")
                f.write("using namespace std;\n")
            for section in sections:
                f.write('\n// %s section\n' % section)
                for slug, desc, ctype, default in PARAMS [section]:
                    if is_extern:
                        f.write('   extern %s %s;\n' %(ctype, slug))
                    else:
                        f.write('%s %s;\n' %(ctype, slug))
           
            if is_extern:
                f.write("}\n")

    def write_pxi(input_filename, output_filename, substitutions):
        data = open(input_filename, 'r').read()
        indent = '    '
        for key, classname, sections in substitutions:
            rep = ''
            for section in sections:
                rep += '%s# %s section\n' % (indent, section)
                for slug, desc, ctype, default in PARAMS[section]:
                    rep += '%sproperty %s:\n' % (indent, slug)
                    rep += '%s    def __get__(%s self):\n' % (indent, classname)
                    rep += '%s        self._check_alive()\n' % indent
                    rep += '%s        return self._thisptr.%s\n' % (indent, slug)
                    rep += '%s\n' % indent
                    rep += '%s    def __set__(%s self, value):\n' % (indent, classname)
                    rep += '%s        self._check_alive()\n' % indent
                    rep += '%s        self._thisptr.%s = <%s> value\n' % (indent, slug, ctype)
                    rep += '%s\n' % indent
            data = data.replace('{{ %s }}' % key, rep)
        with open(output_filename, 'w') as f:
            f.write("""# pxi files\n""")
            f.write("# THIS FILE IS AUTOMATICALLY GENERATED BY param.py\n")
            f.write("# ANY CHANGES WILL BE OVERWRITTEN AT COMPILE TIME\n")
            f.write(data)

    def write_pxd(filename, model, sections):
        with open(filename, 'w') as f:
            # Sadly, Cython docs are incorrect on usage of 'include', so we must include a whole lot of boilerplate
            f.write("""# pxd files\n""")
            f.write("""# THIS FILE IS AUTOMATICALLY GENERATED BY param.py
# ANY CHANGES WILL BE OVERWRITTEN AT COMPILE TIME
from libcpp.vector cimport vector
cdef extern from "c%sModel.h":
    cppclass c%sModel:
        c%sModel() except +  # NB! std::bad_alloc will be converted to MemoryError
        void LoadConfig(const char *filename)
        void Run()
""" % (model, model, model))
            for section in sections:
                f.write('        # %s section\n' % section)
                for slug, desc, ctype, default in PARAMS [section]:
                    f.write('        %s %s\n' % (ctype, slug))

    def write_config_xml(filename, sections):
        with open(filename, 'w') as f:
            f.write("""<!-- THIS FILE IS GENERATED AUTOMATICALLY BY param.py. DO NOT EDIT. -->""")
            f.write('<aturan>')
            for section in sections:
                f.write('    <%s>\n' % section)
                for slug, desc, ctype, default in PARAMS [section]:
                    if ctype == 'bool':
                        default = str(default).lower()  # Python uses True/False, C++, uses true/false
                    f.write('        <%s>%s</%s>  <!-- %s (%s) -->\n' % (slug, default, slug, desc, ctype))
                f.write('    </%s>\n' % section)
            f.write('</aturan>')

    uranus_atmosphere_sections = ['common', 'uranus']

    for filename, classname, sections in [
       ('planet/cUranusDefaults.cpp.inc', 'cUranusModel', uranus_atmosphere_sections),
   ]:
        write_cpp_defaults(filename, classname, sections)

    for filename, classname, sections in [
        ('planet/UranusLoadConfig.cpp.inc', 'cUranusModel', uranus_atmosphere_sections),
   ]:
        write_cpp_load_config(filename, classname, sections)

    for filename, sections in [
        ('planet/UranusParams.h.inc', uranus_atmosphere_sections),
   ]:
        write_cpp_ueaders(filename, sections)

    write_pxi ('python/pyaturan.pyx.template', 'python/pyaturan.pyx', [
        ('uranus_params', 'Uranus', uranus_atmosphere_sections)]
    )

    for filename, model, sections in [
        ('python/uranus_pxd.pxi', 'Uranus', uranus_atmosphere_sections),
   ]:
        write_pxd(filename, model, sections)

    for  filename, sections in [
        ('python/config_aturan.xml', uranus_atmosphere_sections),
   ]:
        write_config_xml(filename, sections)

    for  filename, sections in [
        ('cli/config_aturan.xml', uranus_atmosphere_sections),
   ]:
        write_config_xml(filename, sections)

    for  filename, sections in [
        ('uranus/config_aturan.xml', uranus_atmosphere_sections),
   ]:
        write_config_xml(filename, sections)

if __name__ == '__main__':
    main()
