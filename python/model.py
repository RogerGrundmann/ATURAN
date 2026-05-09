#!/usr/bin/env python

from pyaturan import Uranus

class Model(object):
    """
    ATURAN Uranus' cell structure in the atmosphere
    """
    def __init__(self, cfg_xml):
        self.uran = Uranus()
        self.config_xml = cfg_xml

    def print_config_uran(self, config_xml):
        print("\n\n\n   Uranus' cell structure in the atmosphere") 
        print("\n\n\n   uranus atmosphere configuration file name = ", self.config_xml)
        print("   output path is           ", self.uran.output_path.decode('utf-8'))

    def run_Model_uran(self):
        print("\n   run_Model for the Uranus-Atmosphere code prepared")
        self.uran.run()
        print("\n    successfully terminated Uranus-Atmosphere code")
        print("\n")

uran = Model("config_aturan.xml")
uran.print_config_uran("config_aturan.xml")
uran.run_Model_uran()
