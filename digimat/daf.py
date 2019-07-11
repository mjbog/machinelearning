import os
import numpy as np
from subprocess import PIPE, Popen
import time
import sys
import shutil
import xml.etree.ElementTree as ET
from math import pi
import re
import numpy as np
import operator
#import tensormath as tm
from math import pi, sin, cos
import difflib

try:
    import abaqus
except:
    from scipy.optimize import minimize
    import matplotlib.pyplot as plt
    with_scipy = True
    plotting = True
else:
    with_scipy = False
    plotting = False


from .parsefile import DigimatFile

__author__ = "Michael Bogdanor"

# Load in the digimat MF commands
dirpath = os.path.dirname(__file__)
#print dirpath
MFcommand = dirpath + "\\DigimatMF_cmd.bat"
FABcommand = dirpath + "\\DigimatMF_fabric.bat"
Hybrid_Command = dirpath + "\\launch_convertForHybrid.bat"

idxs = [ '11', '22', '33', '12', '13', '23' ]

############################################################
#  Class for daf files
############################################################
class Daf(DigimatFile):
    """ Object for holding the options from a daf file, initialize the analysis with
      Daf = daf('your-daf-file-name.) """

    def __init__(self,job_name=None,mat_name=False,continuous=False,replace=False):
        """ initialize the dafFile object """

        if job_name is None:
            if continuous:
                job_name = "Continuous"
                shutil.copy(dirpath + "/Continuous_Baseline.mat",job_name+'.daf')
            else:
                job_name = "SFRP_Baseline"
                shutil.copy(dirpath + "/SFRP_Baseline.daf",job_name)
        else:
            if not replace and \
                    (os.path.exists(job_name+'.mat') or os.path.exists(job_name+'.daf')):
                pass
            else:
                if mat_name:
                    shutil.copy(mat_name, job_name+'.mat')
                elif continuous:
                    shutil.copy(dirpath + "/Continuous_Baseline.mat",job_name+'.mat')
                else:
                    shutil.copy(dirpath + "/SFRP_Baseline.daf",job_name + '.mat')

        try:
            if job_name [-4] is '.':
                job_name = job_name[:-4]
        except IndexError:
            pass
        self.job_name = job_name
        self.daf_name = job_name + '.daf'
        self.mat_name = job_name + '.mat'
        self.eng_name = job_name + '.eng'
        self.section_names = []
        self.options = []
        self.values = []
        self.read_daf(mat_name)

    def rename(self,job_name):
        """ Rename the job, .daf, .mat, and .eng files to the given job_name """
        try:
            if job_name [-4] is '.':
                job_name = job_name[:-4]
        except IndexError:
            pass
        self.job_name = job_name
        self.daf_name = job_name + '.daf'
        self.mat_name = job_name + '.mat'
        self.eng_name = job_name + '.eng'

        try:
            file_ext = self.file_name[-4:]
            self.file_name = job_name + file_ext
        except:
            self.file_name = self.mat_name

    ############################################################
    #  Parse the daf file for the options of the analysis
    ############################################################
    def read_daf(self,mat_name=False,savefile=True):
        """ Read in the .daf file and store the information as:
              section_names
               |-> options
               |-> values
        where section_names is an ordered list of the KEYWORD objects
        and options and values are ordered lists of lists corresponding
        to the order of the KEYWORDS in the .daf file with the options and
        values in the order they are written in the .daf file """

        # If mat == False, open the .daf or .mat matching the job_name,
        # otherwise read in the given matfile
        if mat_name:
            if mat_name[-4:] in ['.mat','.daf']:
                file_name = mat_name
            else:
                file_name = mat_name + '.daf'
        else:
            file_name = self.daf_name

        # Try to read in the file as a daf file, otherwise change extension to mat
        try:
            lines = open(file_name).readlines()
        except:
            file_name = file_name.replace('.daf','.mat')
            lines = open(file_name).readlines()

        self.initial_lines = lines
        if savefile:
            self.file_name = file_name
        self.section_names, self.options, self.values = self.read_digimat_file(file_name)

        # Define the isotropy of the fiber
        self.fiber_type = self.get_option('MATERIAL','elastic_model',name='Fiber')

        try:
            self.fabric_name = self.get_option('woven','fabric_file',verbose=False)
        except:
            self.fabric_name = None

        self.read_hybrid_parameters(file_name)

    def read_hybrid_parameters(self, file_name = None):
        """ Read in the computed hybrid parameter section from the given file or the
        matfile for the analysis """
        if file_name == None:
            file_name = self.mat_name

        lines = open(file_name).read()

        try:
            hybrid_param = re.search('# HYBRID(.*\n.*)*',lines).group(0)
        except AttributeError:
            hybrid_param = None

        self.hybrid_param = hybrid_param


    def print_daf(self,new_name = None,verbose=False):
        """shortcut to print_file(self,new_name = None, verbose=False)

        Recreate the digimat file object for the analysis, if no new_name is specified,
        the function will write to the given self.file_name
        verbose = True tells the command to print the name of the file written to"""
        try:
            addenda = self.hybrid_param
        except:
            addenda = ''
        self.print_file(new_name, verbose, addenda)

    ############################################################
    #  Modify the Loading Conditions in the file
    ############################################################
    def loading_uniaxial(self,e11 = 0.03, angle=0):
        """ Modify the mechanical loading case to be one-d"""
        n_sec = self.find_section('LOADING','Mechanical')

        options = [ 'name',
                    'type',
                    'load',
                    'initial_strain',
                    'peak_strain',
                    'history',
                    'quasi_static',
                    'theta_load',
                    'phi_load' ]
        values = [ 'Mechanical',
                   'strain',
                   'uniaxial_1',
                   '0.00000000e+000',
                   '{:16.8e}'.format(e11),
                   'monotonic',
                   'on',
                   '9.000000e+001',
                   '{:16.8e}'.format(angle) ]
        self.options[n_sec] = options
        self.values[n_sec] = values

    def loading_biaxial(self,e11 = 0.03, e22=0.03):
        """ Modify the mechanical loading case to be one-d"""
        if e22 == 0.:
            e22 = 1e-8
        n_sec = self.find_section('LOADING','Mechanical')
        options = [ 'name',
                    'type',
                    'load',
                    'initial_strain',
                    'peak_strain',
                    'strain_ratio',
                    'history',
                    'quasi_static']
        values = [ 'Mechanical',
                   'strain',
                   'biaxial1_2',
                   '0.00000000e+000',
                   '{:16.8e}'.format(e11),
                   '{:16.8e}'.format(e22/e11),
                   'monotonic',
                   'on']
        self.options[n_sec] = options
        self.values[n_sec] = values

    def loading_2D(self,e11 = 0.03, e22 = 0.03, e12 = 0.03):
        """ Modify the mechanical loading case to be general 2-d"""
        n_sec = self.find_section('LOADING','Mechanical')

        options = [ 'name',
                    'type',
                    'load',
                    'initial_strain_11',
                    'strain_11',
                    'initial_strain_22',
                    'strain_22',
                    'initial_strain_12',
                    'strain_12',
                    'history',
                    'quasi_static']
        values = [ 'Mechanical',
                   'strain',
                   'General_2D',
                   '0.00000000e+000',
                   '{:16.8e}'.format(e11),
                   '0.00000000e+000',
                   '{:16.8e}'.format(e22),
                   '0.00000000e+000',
                   '{:16.8e}'.format(e12),
                   'monotonic',
                   'on' ]
        self.options[n_sec] = options
        self.values[n_sec] = values


    ############################################################
    #  Run and Analyze the results
    ############################################################
    def run_MF(self,verbose=0):
        """ Run the given Digimat .mat file """
        self.print_daf(self.mat_name)
        matAbs = os.path.abspath(self.mat_name)
        cmd = MFcommand+' input="'+matAbs+'"'
        if verbose >= 1:
            print(cmd)
        subp = Popen([MFcommand,'input='+matAbs], stdout=PIPE, stderr=PIPE)
        out, err = subp.communicate()

        if verbose >= 2:
            print(out)

        if verbose >= 3:
            print(err)

        return str(out),str(err)

    def read_eng(self,verbose=0):
        """ Read the .eng file from the analysis to retrieve homogenized
        material properties """
        try:
            # print 'Read info from', self.eng_name
            lines = open(self.eng_name).read()
        except (IOError):
            print("Could not find file named", self.eng_name)
            return None

        results = {}

        # Compliance
        C = _strip_data(re.search('Compliance.*\n(.*\n)+?#',lines).group())
        results['C'] = C

        Ex = 1/C[0,0]
        Ey = 1/C[1,1]
        Ez = 1/C[2,2]
        Gxy = 1/(C[3,3])
        Gyz = 1/(C[4,4])
        Gxz = 1/(C[5,5])
        vyx = -Ey*C[0,1]
        vxy = -Ex*C[1,0]
        vzx = -Ez*C[0,2]
        vxz = -Ex*C[2,0]
        vzy = -Ez*C[1,2]
        vyz = -Ey*C[2,1]
        if verbose > 0:
            print('E1:  {0:6.0f}   E2:  {1:6.0f}   E3:  {2:6.0f}'.format(Ex, Ey, Ez))
            print('G12: {0:6.0f}   G13: {1:6.0f}   G23: {2:6.0f}'.format(Gxy, Gxz, Gyz))
            print('v12: {0:6.4f}   v13: {1:6.4f}   v23: {2:6.4f}'.format(vxy, vxz, vyz))

        results['Ex'] = Ex
        results['Ey'] = Ey
        results['Ez'] = Ez
        results['Gxy'] = Gxy
        results['Gxz'] = Gxz
        results['Gyz'] = Gyz
        results['vxy'] = vxy
        results['vxz'] = vxz
        results['vyz'] = vyz

        # Stiffness
        K = _strip_data(re.search('Stiffness.*\n(.*\n)+?#',lines).group())
        results['K'] = K

        # CTE
        try:
            A = _strip_data(re.search('Thermal Expansion.*\n(.*\n)+?#',lines).group())[0]
        except AttributeError:
            pass
        else:
            alpha = np.array( [ [ A[0], A[5], A[4] ],
                                [ A[5], A[1], A[3] ],
                                [ A[4], A[3], A[2] ] ])
            results['cte'] = alpha

        # Density
        try:
            density = float(re.search('Global density.*\n.*\n.*',lines).group().splitlines()[-1])
        except AttributeError:
            density = None
        else:
            results['density'] = density
            if verbose:
                print('Density: {0:6.4f}'.format(density))

        self.results = results
        self.density = density
        self.Es = [ Ex, Ey, Ez]
        self.Gs = [ Gxy, Gxz, Gyz ]
        self.vs = [ vxy, vxz, vyz ]
        self.cte = alpha

        return results

    def read_mac(self,idx=1):
        """ Read the macroscopic stress-strain data for failure, returns the stress
        for the given loading index"""
        mac = self.job_name + '.mac'
        lines = open(mac).readlines()

        time = []
        s11 = []
        s = []
        e11 = []
        e = []
        temp = []

        for line in lines:
            try:
                data = [ float(x) for x in line.split()]
                time.append(data[0])
                e.append(data[1:7])
                s.append(data[8:14])
                e11.append(data[idx])
                s11.append(data[idx+7])
                if len(data) > 15:
                    temp.append(data[15])
            except:
                pass
        self.s11 = s11
        self.e11 = e11
        self.e = e
        self.s = s
        self.time = time
        self.temp = temp
        self.maxs = max(s11)
        return s11, e11

    def run_analysis(self,loglevel = 0):
        """ Run the mean field homogenization analysis and return the results from
        the .eng and .mac files """
        self.print_daf(self.daf_name)
        self.print_daf(self.mat_name)
        out,err = self.run_MF(loglevel)

        while 'License check out failed' in out:
            print('DIGIMAT >>> License check out failed, wait 5 sec')
            time.sleep(5)
            out,err = self.run_MF(loglevel)
            print(out)
        if 'Unable to open file' in out:
            print('DIGIMAT >>> File could not be accessed, wait 5 sec')
            time.sleep(5)
            out,err = self.run_MF(loglevel)
        elif 'ERROR' in out:
            err_lines = err.split('\n')
            for line in err_lines:
                if 'ERROR' in line:
                    print('DIGIMAT >>>', line)
            out_lines = out.split('\n')
            for line in out_lines:
                if 'ERROR' in line:
                    print('DIGIMAT >>>', line)
            return None

        results = self.read_eng()
        #self.read_mac()
        return results

    def read_dfe(self):
        """ Read through the failure envelope file """
        dfe = self.job_name + '_FailureEnvelope.dfe'
        lines = open(dfe).readlines()
        FE = []
        for line in lines:
            try:
                data = [ float(x) for x in line.split() ]
                FE.append(data)
            except:
                pass
        self.FE = FE
        return FE

    def set_fiber_type(self,fiber_type):
        """ Change the isotropy of the fiber, fiber type can be either 'trans(versely isotropic)' or 'isotropic'. If
        fiber type contains 'trans', the fiber will be set to transversely isotropic, else it will be set to isotropic"""
        # Find the fiber options
        current_model = self.get_option('MATERIAL','elastic_model',name='Fiber')

        if 'trans' in fiber_type:
            new_model = 'transversely_isotropic'
            self.fiber_type = 'transversely_isotropic'
        else:
            new_model = 'isotropic'
            self.fiber_type = 'isotropic'

        self.set_option('MATERIAL','elastic_model',new_model,'Fiber')

        sec_no = None
        # Identify the fiber material section
        for i in range(len(self.values)):
            if self.section_names[i] == 'MATERIAL' and self.get_option(i,'name') == 'Fiber':
                sec_no = i

        if sec_no is None:
            print('No section for MATERIAL named Fiber could be found.')

        options = self.options[sec_no]
        values = self.values[sec_no]

        if new_model == 'isotropic':
            for i in range(len(options))[::-1]:
                # Delete inPlane_Young
                if options[i] == 'inPlane_Young':
                    del options[i]
                    del values[i]
                # Delete transverse_shear
                elif options[i] == 'transverse_shear':
                    del options[i]
                    del values[i]
                # Change axial Young to Young
                elif options[i] == 'axial_Young':
                    if 'Young' in options:
                        del options[i]
                        del values[i]
                    else:
                        options[i] = 'Young'
                # Delete transverse Poisson
                elif options[i] == 'transverse_Poisson':
                    del options[i]
                    del values[i]
                # Change inPlane Poisson to Poisson
                elif options[i] == 'inPlane_Poisson':
                    if 'Poisson' in options:
                        del options[i]
                        del values[i]
                    else:
                        options[i] = 'Poisson'
                # Delete in plane CTE
                elif options[i] == 'inPlane_CTE':
                    del options[i]
                    del values[i]
                # Change axial CTE to thermal_expansion
                elif options[i] == 'axial_CTE':
                    options[i] = 'thermal_expansion'

        else:
            # Find the youngs modulus and add in the sections for inplane and transverse moduli
            for i in range(len(options)):
                if options[i] == 'Young':
                    E1 = float(values[i])
                    options[i] = 'axial_Young'
                    options.insert(i,'inPlane_Young')
                    options.insert(i+1,'transverse_shear')
                    values.insert(i,'{0:16.8e}'.format(E1*0.1))
                    values.insert(i+1,'{0:16.8e}'.format(E1*0.15))
                    break
            # Add in the transverse Poisson
            for i in range(len(options)):
                if options[i] == 'Poisson':
                    options[i] = 'inPlane_Poisson'
                    options.insert(i,'transverse_Poisson')
                    values.insert(i,values[i])
                    break
            # Add in plane CTE
            for i in range(len(options)):
                if options[i] == 'thermal_expansion':
                    options[i] = 'axial_CTE'
                    options.insert(i,'inPlane_CTE')
                    values.insert(i,values[i])
                    break

        self.options[sec_no] = options
        self.values[sec_no] = values
        self.print_daf()


    # Important options:
    # FIBER
    #   MATERIAL, name = Fiber
    #     density
    #     Young
    #     Poisson
    #     thermal_expansion
    #   PHASE, name = Fiber
    #     mass_fraction
    #     orientation_11
    #     orientation_22
    #     orientation_33
    #     orientation_12
    #     orientation_13
    #     orientation_23
    #   FAILURE INDICATOR, name = Fiber
    #     tensile_strength
    #     compressive_strength
    def _set_get_fiber(self, set_=True, **kwargs):
        """ Define the fiber properties based on the following kwargs:
          density      kg/m^3
          modulus      MPa
          poisson      unitless
          cte          strain/degC
          mass_frac
          vol_frac
          orientation  scalar, diagonal, or full 2D tensor
          t_strength    MPa
          c_strength    MPa
          aspect_ratio
        """
        if set_ == True:
            func = self.set_option
        else:
            func = self.get_option

        if kwargs is not None:
            for key, value in kwargs.items():
                if not set_:
                    value = 0
                if key.lower() == 'density':
                    value = func('MATERIAL','density','{:16.8e}'.format(value),'Fiber')
                elif key.lower() == 'modulus':
                    if 'transverse' in self.get_option('MATERIAL','elastic_model',name='Fiber'):
                        try:
                            E1 = float(value)
                            E2 = float(value)*0.1
                            G12 = float(value)*0.15
                        except TypeError:
                            E1, E2, G12 = value
                        E2 = func('MATERIAL','inPlane_Young','{:16.8e}'.format(E2),'Fiber')
                        G12 = func('MATERIAL','transverse_shear','{:16.8e}'.format(G12),'Fiber')
                        E1 = func('MATERIAL','axial_Young','{:16.8e}'.format(E1),'Fiber')
                        value = [E1, E2, G12]
                    else:
                        value = func('MATERIAL','Young','{:16.8e}'.format(value),'Fiber')
                elif key.lower() == 'poisson':
                    if 'transverse' in self.get_option('MATERIAL','elastic_model',name='Fiber'):
                        try:
                            v12 = float(value)
                            v23 = float(value)
                        except TypeError:
                            v12, v23 = value
                        v12 = func('MATERIAL','transverse_Poisson','{:16.8e}'.format(v12),'Fiber')
                        v23 = func('MATERIAL','inPlane_Poisson','{:16.8e}'.format(v23),'Fiber')
                        value = [v12, v23]
                    else:
                        value = func('MATERIAL','Poisson','{:16.8e}'.format(value),'Fiber')
                elif key.lower() == 'cte':
                    if 'transverse' in self.get_option('MATERIAL','elastic_model',name='Fiber'):
                        try:
                            cte1 = float(value)
                            cte2 = float(value)
                        except TypeError:
                            cte1, cte2 = value
                        cte1 = func('MATERIAL','axial_CTE','{:16.8e}'.format(cte1),'Fiber')
                        cte2 = func('MATERIAL','inPlane_CTE','{:16.8e}'.format(cte2),'Fiber')
                        value = [cte1, cte2]
                    else:
                        value = func('MATERIAL','thermal_expansion','{:16.8e}'.format(value),'Fiber')
                elif key.lower() == 'mass_frac':
                    if set_:
                        self.delete_option('PHASE','volume_fraction','Fiber')
                        self.delete_option('PHASE','volume_fraction','Matrix')
                        func('PHASE','mass_fraction','{:16.8e}'.format(1.-value),'Matrix')
                    value = func('PHASE','mass_fraction','{:16.8e}'.format(value),'Fiber')
                elif key.lower() == 'vol_frac':
                    if set_:
                        self.delete_option('PHASE','mass_fraction','Fiber')
                        self.delete_option('PHASE','mass_fraction','Matrix')
                        try:
                            v_vol = float(self.get_option('PHASE','volume_fraction',0,'Voids'))
                        except:
                            v_vol = 0.
                        func('PHASE','volume_fraction','{:16.8e}'.format(1.-value-v_vol),'Matrix')
                    value = func('PHASE','volume_fraction','{:16.8e}'.format(value),'Fiber')
                elif key.lower().startswith('orient'):
                    if set_:
                        self.set_fiber_orientation(value)
                    else:
                        value = self.get_fiber_orientation()
                elif key.lower() == 't_strength':
                    value = func('FAILURE INDICATOR','tensile_strength','{:16.8e}'.format(value),'Fiber')
                elif key.lower() == 'c_strength':
                    value = func('FAILURE INDICATOR','compressive_strength','{:16.8e}'.format(value),'Fiber')
                elif 'aspect' in key.lower():
                    value = func('PHASE','aspect_ratio','{:5f}'.format(value),'Fiber')
                else:
                    print('Could not recognize key:', key)

        if not set_:
            return value

    def set_fiber(self, **kwargs):
        """ Define the fiber properties based on the following kwargs:
          density      kg/m^3
          modulus      MPa
          poisson      unitless
          cte          strain/degC
          mass_frac
          vol_frac
          orientation  scalar, diagonal, or full 2D tensor
          t_strength    MPa
          c_strength    MPa
          aspect_ratio
        """
        self._set_get_fiber(True, **kwargs)

    def get_fiber(self, *args):
        if args is not None:
            kwargs = {}
            values = {}
            for key in args:
                kwargs[key] = None
                value = self._set_get_fiber(False, **kwargs)
                try:
                    value = float(value)
                except:
                    pass
                values[key] = value
            if len(values) > 1:
                return values
            else:
                return list(values.values())[0]
            return values

    def get_fiber_orientation(self):
        """ Print the orientation tensor for the current analysis or the orientation file"""
        self.orient_type = self.get_option('PHASE','orientation',name='Fiber')
        if self.orient_type == 'tensor':
            idxs = [ '11', '22', '33', '12', '13', '23' ]
            tensor = np.zeros((3,3))
            for i in range(6):
                value = float(self.get_option('PHASE','orientation_'+idxs[i],name='Fiber'))
                k = int(idxs[i][0])-1
                l = int(idxs[i][1])-1
                tensor[k,l] = value
                tensor[l,k] = value
            self.orientation = tensor
            return tensor
        else:
            self.orientation_file = self.get_option('PHASE','orientation_file',name='Fiber')
            return self.orientation_file

    def set_fiber_orientation(self, orientation):
        """ edit the phase orientation options to evalueuate a different orientation
        state in the fiber, can be passed in as a:
          3x3 matrix - full 2nd order orientation tensor
          3x1 vector - diagonal of the orientation tensor
          scalar - first component of the in-plane orientation tensor
        """
        #print 'test'
        if type(orientation) == type('str'):
            #print 'process as string'
            for i in range(6):
                self.delete_option('PHASE','orientation_'+idxs[i],name='Fiber')
            self.set_option('PHASE','orientation','file',name='Fiber')
            self.set_option('PHASE','orientation_format','DIGIMAT',name='Fiber')
            self.set_option('PHASE','orientation_file',orientation,name='Fiber')

            # Set the analysis options to use Abaqus instead of Digimat
            self.set_option('ANALYSIS', '')

        else:
            #print 'process as number'
            self.set_option('PHASE','orientation','tensor',name='Fiber')
            self.delete_option('PHASE','orientation_format',name='Fiber')
            self.delete_option('PHASE','orientation_file',name='Fiber')

            # Set the analysis options to use Digimat instead of Abaqus

            orient_tensor = np.zeros(6)
            try:
                orient_tensor = [ orientation[0][0],
                                  orientation[1][1],
                                  orientation[2][2],
                                  orientation[0][1],
                                  orientation[0][2],
                                  orientation[1][0] ]
            except:
                try:
                    for i in range(len(orientation)):
                        orient_tensor[i] = orientation[i]
                except:
                    orient_tensor = [ orientation, 1.0-orientation, 0., 0., 0., 0. ]
            for i in range(6):
                self.set_option('PHASE','orientation_'+idxs[i],
                                '{:12.4e}'.format(orient_tensor[i]),'Fiber')

    def set_digimat_loading(self):
        """ Create a section for the digimat thermal and mechanical loading if not already defined """




    def remove_digimat_loading(self):
        """ Delete the loading sections from the analysis"""
        self.delete_section('LOADING', 'Mechanical')
        self.delete_section('LOADING', 'Thermal')




    def get_fiber_volume_fraction(self):
        """ return the fiber volume fraction defined in the material or from the density and mass fraction """
        vol_frac = self.get_fiber('vol_frac')
        if vol_frac is None:
            mass_frac = self.get_fiber('mass_frac')
            d_fiber = self.get_fiber('density')
            d_matrix = self.get_matrix('density')
            vol_frac = ( mass_frac*d_matrix ) / ( d_fiber + mass_frac*( d_matrix - d_fiber))
        return vol_frac

    def get_fiber_mass_fraction(self):
        """ return the fiber mass fraction defined in the material or from the density and volume fraction """
        mass_frac = self.get_fiber('mass_frac')
        vol_frac = self.get_fiber('vol_frac')
        d_fiber = self.get_fiber('density')
        d_matrix = self.get_matrix('density')
        if mass_frac is None:
            m_fiber = vol_frac*d_fiber
            m_matrix = (1. - vol_frac)*d_matrix
            mass_frac = m_fiber/(m_fiber + m_matrix)
        return mass_frac

    def get_composite_density(self):
        """ Compute the density of the composite """
        try:
            d_fiber = self.get_fiber('density')
        except:
            print('No fiber density given')
        try:
            d_matrix = self.get_matrix('density')
        except:
            print('No matrix density given')
        vol_frac = self.get_fiber('vol_frac')
        if vol_frac is None:
            vol_frac = self.get_fiber_volume_fraction()
        self.composite_density = vol_frac*d_fiber + (1-vol_frac)*d_matrix
        return self.composite_density

    # MATRIX
    #   MATERIAL, name = Matrix
    #     density
    #     Young
    #     Poisson
    #     thermal_expansion
    #   PHASE, name = Matrix
    #     mass_fraction
    #   FAILURE INDICATOR, name = Matrix
    #     tensile_strength
    #     compressive_strength
    def _set_get_matrix(self, set_=True, **kwargs):
        """ Define the matrix properties based on the following kwargs:
          density      kg/m^3
          modulus      MPa
          poisson      unitless
          cte          strain/degC
          mass_frac - set using the defineFiber method
          vol_frac - set using the defineFiber method
          t_strength    MPa
          c_strength    MPa
        """

        if set_ == True:
            func = self.set_option
        else:
            func = self.get_option

        if kwargs is not None:
            for key, value in kwargs.items():
                if not set_:
                    value = 0

                if key.lower() == 'density':
                    value = func('MATERIAL','density','{:16.8e}'.format(value),'Matrix')
                elif key.lower() == 'modulus':
                    value = func('MATERIAL','Young','{:16.8e}'.format(value),'Matrix')
                elif key.lower() == 'poisson':
                    value = func('MATERIAL','Poisson','{:16.8e}'.format(value),'Matrix')
                elif key.lower() == 'cte':
                    value = func('MATERIAL','thermal_expansion','{:16.8e}'.format(value),'Matrix')
                elif key.lower() == 't_strength':
                    value = func('FAILURE INDICATOR','tensile_strength','{:16.8e}'.format(value),'Matrix')
                elif key.lower() == 'c_strength':
                    value = func('FAILURE INDICATOR','compressive_strength','{:16.8e}'.format(value),'Matrix')
        if not set_:
            return value

    def set_matrix(self, **kwargs):
        """ Define the matrix properties based on the following kwargs:
          density      kg/m^3
          modulus      MPa
          poisson      unitless
          cte          strain/degC
          mass_frac - set using the defineFiber method
          vol_frac - set using the defineFiber method
          t_strength    MPa
          c_strength    MPa
        """
        self._set_get_matrix(True, **kwargs)

    def get_matrix(self, *args):
        if args is not None:
            kwargs = {}
            values = {}
            for key in args:
                kwargs[key] = None
                value = self._set_get_matrix(False, **kwargs)
                try:
                    value = float(value)
                except:
                    pass
                values[key] = value
            if len(values) > 1:
                return values
            else:
                return list(values.values())[0]

    def set_matrix_failure(self,failure_type,tensile_strength,
                           compressive_strength=None,axial_inplane=1e-8,
                           shear_strength=None, verbose=0):
        """ Update the failure indicator properties for the matrix for the given
        failure type and failure values """
        # For the given failure type, set the appropriate parameters in the failure
        # indicator
        matrix_variable_names = [ 'axial_tensile_strength',
                                  'axial_compressive_strength',
                                  'inplane_tensile_strength',
                                  'inplane_compressive_strength',
                                  'transverse_shear_strength',
                                  'axial_inplane_strength',
                                  'inplane_shear_strength',
                                  'tensile_strength',
                                  'compressive_strength' ]

        matrix_variable_values = { 'axial_tensile_strength':1,
                                   'axial_compressive_strength':2,
                                   'inplane_tensile_strength':1,
                                   'inplane_compressive_strength':2,
                                   'transverse_shear_strength':3,
                                   'axial_inplane_strength':0,
                                   'inplane_shear_strength':3,
                                   'tensile_strength':1,
                                   'compressive_strength':2 }

        failure_indicators = { 'Tsai_Hill_2D'      :[  1, 1, 1, 1, 1, 0, 0, 0, 0],
                               'Azzi_Tsai_Hill_2D' :[  1, 1, 1, 1, 1, 0, 0, 0, 0],
                               'Tsai_Wu_2D'        :[  1, 1, 1, 1, 1, 1, 0, 0, 0],
                               'maximum_stress_2D' :[  1, 1, 1, 1, 1, 0, 0, 0, 0],
                               'Hashin_Rotem_2D'   :[  1, 1, 1, 1, 1, 0, 0, 0, 0],
                               'Hashin_2D'         :[  1, 1, 1, 1, 1, 0, 1, 0, 0],
                               'Hashin_3D'         :[  1, 1, 1, 1, 1, 0, 1, 0, 0],
                               'Christensen'       :[  0, 0, 0, 0, 0, 0, 0, 1, 1] }

        # Check that the failure type is in the failure indicators
        if not failure_type in failure_indicators.keys():
            if verbose > 0:
                print('Could not find failure type: ', failure_type)
            close_matches = difflib.get_close_matches(failure_type,failure_indicators.keys())
            if len(close_matches) > 0:
                if verbose > 0:
                    print('Using', close_matches[0], 'instead')
                failure_type = close_matches[0]
            else:
                print('No closely matching failure types found')
                raise

        # Basic options
        new_options = []
        new_values = []
        new_options.append('name')
        new_values.append('Matrix')
        new_options.append('type')
        new_values.append(failure_type)
        new_options.append('use_linear_formulation')
        new_values.append('on')
        new_options.append('axes')
        new_values.append('local')

        if compressive_strength == None:
            compressive_strength = tensile_strength
        if shear_strength == None:
            shear_strength = tensile_strength

        failure_values = [ axial_inplane,
                           tensile_strength,
                           compressive_strength,
                           shear_strength ]

        for idx, fi_flag in enumerate(failure_indicators[failure_type]):
            if fi_flag:
                variable = matrix_variable_names[idx]
                new_options.append(variable)
                new_values.append(failure_values[matrix_variable_values[variable]])

        section_i = self.find_section('FAILURE INDICATOR','Matrix')
        self.options[section_i] = new_options
        self.values[section_i] = new_values

    ### TODO Additional Matrix Considerations

    def convert_to_strain_failure(self, nonlinearity=1.):
        """ Convert the failure strengths in the material model to failure strains
        using the elastic properties of the material. Nonlinearity increases the
        failure strain beyond the elastic behavior """
        # Convert the fiber
        fib_sec = self.find_section('Failure Indicator','Fiber')
        opts = self.options[fib_sec]
        vals = self.values[fib_sec]

        E_fiber = float(self.get_option('material','axial_Young',name='Fiber'))

        type_idx = opts.index('type')
        vals[type_idx] = 'maximum_strain'
        X1T = float(vals[opts.index('tensile_strength')])
        if X1T > 1:
            vals[opts.index('tensile_strength')] = '{0:5.4e}'.format(X1T/E_fiber)

        X1C = float(vals[opts.index('compressive_strength')])
        if X1C > 1:
            vals[opts.index('compressive_strength')] ='{0:5.4e}'.format(X1C/E_fiber)

        # Convert the Matrix
        mat_sec = self.find_section('Failure Indicator','Matrix')
        opts = self.options[mat_sec]
        vals = self.values[mat_sec]

        E_matrix = float(self.get_option('material','Young',name='Matrix'))
        v_matrix = float(self.get_option('material','Poisson',name='Matrix'))
        G_matrix = E_matrix/ ( 2.*(1. + v_matrix))

        type_idx = opts.index('type')
        vals[type_idx] = 'Tsai_Wu_3D_transversely_isotropic_Strain'

        names = [ '{0:1s}_{1:1s}_strength'.format(a,b) for a in ['axial','inplane'] for b in ['tensile','compressive'] ]

        for name in names:
            X = float(vals[opts.index(name)])
            if X > 1:
                vals[opts.index(name)] = '{0:5.4e}'.format(X/E_matrix)

        X = float(vals[opts.index('transverse_shear_strength')])
        if X > 1:
            vals[opts.index('transverse_shear_strength')] = '{0:5.4e}'.format(X/G_matrix)

        try:
            self.delete_option('Failure Indicator','inplane_shear_strength','Matrix')
        except:
            pass



    ### TODO Progressive Damage









    ############################################################
    #  Voids
    ############################################################
    def add_voids(self,vol_frac,aspect_ratio = 10, orientation=1.0):
        """ Add a phase to the analysis containing voids """
        # Define the structure of the orientation tensor
        orient_tensor = np.zeros(6)
        try:
            orient_tensor = [ orientation[0][0],
                              orientation[1][1],
                              orientation[2][2],
                              orientation[0][1],
                              orientation[0][2],
                              orientation[1][0] ]
        except:
            try:
                for i in range(len(orientation)):
                    orient_tensor[i] = orientation[i]
            except:
                orient_tensor = [ orientation, 1.0-orientation, 0., 0., 0., 0. ]

        # Create the structure of the void phase
        v_options = [ 'name',
                      'type',
                      'volume_fraction',
                      'behavior',
                      'aspect_ratio',
                      'orientation',
                      'orientation_11',
                      'orientation_22',
                      'orientation_33',
                      'orientation_12',
                      'orientation_13',
                      'orientation_23',
                      'closure',
                      'coated' ]
        v_values = [ 'Voids',
                     'void',
                     '{:15.3e}'.format(max(vol_frac,1e-6)),
                     'deformable_solid',
                     '{:15.3e}'.format(aspect_ratio),
                     'tensor'] + \
                   ['{:15.3e}'.format(o) for o in orient_tensor] + \
                   ['orthotropic',
                    'no']

        # Check if the daf has a phase of type void
        void_sec = None
        last_phase = None
        for idx, name in enumerate(self.section_names):
            if name.lower() == 'phase':
                last_phase = idx
                for opt,val in zip(self.options[idx],self.values[idx]):
                    if opt.lower() == 'type' and val.lower() == 'void':
                        void_sec = idx
        if void_sec == None:
            self.section_names = self.section_names[:last_phase] + ['PHASE'] + \
                                 self.section_names[last_phase:]
            self.options = self.options[:last_phase] + [v_options] + \
                           self.options[last_phase:]
            self.values = self.values[:last_phase] + [v_values] + \
                          self.values[last_phase:]
        else:
            self.options[void_sec] = v_options
            self.values[void_sec] = v_values

        # Check the MICROSTRUCTURE phase
        for idx, name in enumerate(self.section_names):
            if name.lower() == 'microstructure':
                phases = []
                for opt, val in zip(self.options[idx],self.values[idx]):
                    if opt.lower() == 'phase':
                        phases.append(val)

                if 'Voids' in phases:
                    pass
                else:
                    self.options[idx].append('phase')
                    self.values[idx].append('Voids')

        self._verify_volume_fraction()

    def remove_voids(self):
        """ Remove the void phase """
        void_sec = None
        last_phase = None
        for idx, name in enumerate(self.section_names):
            if name.lower() == 'phase':
                last_phase = idx
                for opt,val in zip(self.options[idx],self.values[idx]):
                    if opt.lower() == 'type' and val.lower() == 'void':
                        void_sec = idx
        if void_sec == None:
            pass
        self._verify_volume_fraction()

    def _verify_volume_fraction(self):
        """ After adding or removing voids, make sure the volume fractions in the
        constituents match up """
        f_vf = self.get_fiber_volume_fraction()
        self.set_fiber(vol_frac=f_vf)

    ############################################################
    #  Define different loading conditions
    ############################################################
    def help_loading(self):
        print('set_load_orientation(phi)')
        print(self.set_load_orientation.__doc__)
        print('set_peak_strain(e)')
        print(self.set_peak_strain.__doc__)
        print('set_compression(e_c)')
        print(self.set_compression.__doc__)
        print('set_axial_loading(phi)')
        print(self.set_axial_loading.__doc__)
        print('set_transverse_loading(phi)')
        print(self.set_transverse_loading.__doc__)
        print('set_shear_loading(phi)')
        print(self.set_shear_loading.__doc__)

    def set_load_orientation(self,phi):
        """ Change the angle at which the load is applied to the RVE """
        self.set_option('LOADING','phi_load','{:10.4e}'.format(phi))

    def set_peak_strain(self,e = 0.05):
        """ Change the peak strain in the loading to the given value """
        self.set_option('LOADING','peak_strain','{:16.8e}'.format(e))

    def set_compression(self,e_c = -0.05):
        """ Change the loading case to compression (forces negative value)"""
        self.set_option('LOADING','peak_strain','{:16.8e}'.format(-abs(e_c)))

    def set_axial_loading(self, e = None):
        """ Change the loading case to be applied axially, optional argument for strain """
        self.set_option('LOADING','load','uniaxial_1')
        self.set_option('LOADING','phi_load','0.')
        if not e == None:
            self.set_peak_strain(e)

    def set_transverse_loading(self, e = None):
        """ Change the loading case to be applied transverse to the fiber, optional
        argument for strain  """
        self.set_option('LOADING','load','uniaxial_1')
        self.set_option('LOADING','phi_load','90.')
        if not e == None:
            self.set_peak_strain(e)

    def set_shear_loading(self, e = None):
        """ Change the loading case to be applied in shear, optional argument for
        strain  """
        self.set_option('LOADING','load','shear_12')
        self.set_option('LOADING','phi_load','0.')
        if not e == None:
            self.set_peak_strain(e)

    def set_composite_failure(self, X1T, X1C, X2T, X2C, XIS, XTS):
        """ set_composite_failure(X1T, X1C, X2T, X2C, XIS, XTS)
          Create a composite failure indicator for the material for Hashin 3D with
          progressive failure """

        # Delete the failure indicators from the fiber and matrix
        self.delete_option('MATERIAL','failure_indicator','Fiber')
        self.delete_option('MATERIAL','failure_indicator','Matrix')

        # Reset failure options in the analysis
        self.set_option('ANALYSIS','PF_max_damage','{0:8.6e}'.format(0.99))
        self.set_option('ANALYSIS','PF_critical_damage','none')
        self.set_option('ANALYSIS','PF_explicit_formulation','on')
        self.set_option('ANALYSIS','failure_indicator','Composite')

        # Create the new composite failure indicator

        self.section_names.append('FAILURE INDICATOR')

        opt_vals = [ [ 'name', 'Composite' ],
                     [ 'type', 'Hashin_3D' ],
                     [ 'use_linear_formulation', 'on' ],
                     [ 'axes', 'local' ],
                     [ 'axial_tensile_strength', '{0:16.8e}'.format(X1T) ],
                     [ 'axial_compressive_strength', '{0:16.8e}'.format(X1C) ],
                     [ 'inplane_tensile_strength', '{0:16.8e}'.format(X2T) ],
                     [ 'inplane_compressive_strength', '{0:16.8e}'.format(X2C) ],
                     [ 'transverse_shear_strength', '{0:16.8e}'.format(XTS) ],
                     [ 'inplane_shear_strength', '{0:16.8e}'.format(XIS) ] ]
        opts = [ opt_val[0] for opt_val in opt_vals ]
        vals = [ opt_val[1] for opt_val in opt_vals ]

        self.options.append(opts)
        self.values.append(vals)

        # Create the new progressive damage model
        self.section_names.append('PROGRESSIVE FAILURE')

        opt_vals = [ [ 'name', 'PFM_Composite' ],
                     [ 'type', 'Matzenmiller_3D' ],
                     [ 'failure_indicator', 'Composite' ],
                     [ 'damage_law_fA', 'heaviside,1.00,1.00' ],
                     [ 'damage_law_fB', 'heaviside,1.00,1.00' ],
                     [ 'damage_law_fC', 'exponential,1.00,1.00,5.00,1.00' ],
                     [ 'damage_law_fD', 'exponential,1.00,1.00,5.00,1.00' ] ]
        opts = [ opt_val[0] for opt_val in opt_vals ]
        vals = [ opt_val[1] for opt_val in opt_vals ]

        self.options.append(opts)
        self.values.append(vals)

    ############################################################
    #  Analysis Options
    ############################################################
    # ANALYSIS
    #   homogenization_model
    #     Mori_Tanaka
    #     double_inclusion
    #   number_angle_increments
    def _set_get_analysis_options(self, set_=True, **kwargs):
        """ Define the analysis options:
          mori_tanaka  True/False
          double       True/False
          angle_inc    # of angle increments
        """

        if set_ == True:
            func = self.set_option
        else:
            func = self.get_option

        if kwargs is not None:
            for key, value in kwargs.items():
                if not set_:
                    value = '0'
                if key.lower().startswith('mori'):
                    if value:
                        value = func('ANALYSIS','homogenization_model','Mori_Tanaka')
                    else:
                        value = func('ANALYSIS','homogenization_model','double_inclusion')
                if key.lower().startswith('double'):
                    if value:
                        value = func('ANALYSIS','homogenization_model','double_inclusion')
                    else:
                        value = func('ANALYSIS','homogenization_model','Mori_Tanaka')
                if key.lower() == 'angle_inc':
                    value = func('ANALYSIS','number_angle_increments',str(int(value)))
            return value

    def set_analysis_options(self, **kwargs):
        """ Define the analysis options:
          mori_tanaka  True/False
          double       True/False
          angle_inc    # of angle increments
        """
        self._set_get_analysis_options(True, **kwargs)

    def get_analysis_options(self, *args):
        if args is not None:
            kwargs = {}
            values = {}
            for key in args:
                kwargs[key] = None
                value = self._set_get_analysis_options(False, **kwargs)
                values[key] = value
                try:
                    value = float(value)
                except:
                    pass
            if len(values) > 1:
                return values
            else:
                return values.values()[0]

    def _clean_FE_options(self):
        """ Remove all options in the ANALYSIS section that start with DPP """
        sec_num = self.find_section('ANALYSIS')
        options = self.options[sec_num]
        values = self.values[sec_num]
        for i in reversed(range(len(options))):
            if options[i].startswith('DPP'):
                #print options[i], values[i]
                del options[i]
                del values[i]

    ############################################################
    #  Process the output information from the Digimat Analysis
    ############################################################
    def get_stress_strain(self, plot=False):
        """ Read the mac file and print the stress-strain result """
        s,e = self.read_mac()
        self.s = s
        self.e = e
        idx, self.maxs = max(enumerate(s), key=operator.itemgetter(1))
        self.maxe = e[idx]

        if plot and plotting:
            print('Max stress =', max(s))
            plt.plot(e,s)
            plt.show()

    def get_eng_constants(self):
        """ Read the current .eng file for engineering constants """
        r = self.read_eng()
        constants = [ r['Ex'], r['Ey'], r['e_z'],
                      r['vxy'], r['vxz'], r['vyz'],
                      r['Gxy'], r['Gxz'], r['Gyz'] ]
        return constants

    ############################################################
    #  Calibration Methods
    ############################################################
    def set_matrix_density(self,composite_density):
        """ Calibrate the density of the matrix to the density of the composites based on the known fiber density """
        self.composite_density = composite_density
        vol_frac = self.get_fiber('vol_frac')
        d_fiber = self.get_fiber('density')
        if vol_frac is None:
            # Use volume fraction
            mass_frac = self.get_fiber_mass_fraction()
            d_matrix = (( 1 - mass_frac)*composite_density / \
                        ( 1 - mass_frac*composite_density/d_fiber ))
        else:
            d_matrix = (composite_density - vol_frac*d_fiber)/(1-vol_frac)
        self.set_matrix(density=d_matrix)
        return d_matrix

    ############################################################
    #  Hybrid Parameter Computation
    ############################################################
    def compute_hybrid_parameters(self):
        """ This method computes the hybrid parameters for the given material file
          calls the launch_convertForHybrid.bat command in the python directory

          Arguments to the batch file:
            1. executable
            2. configuration file
            3. testcases location (with .mat file)
            4. number of cpus
            5. hybrid_param location

        """

        if not self.hybrid_param == None:
            print(
                'Hybrid parameters have already been computed for this analyis. To recompute, set self.hybrid_param=None run self.print_daf, and try again')
            return

        # The executable and configuration file are located in the python directory
        exec_name = dirpath + '\\convertForHybrid_v1.2.4\\convertForHybrid_v1.2.4\\convert_for_hybrid.exe'
        config_file = dirpath + '\\configs\\example.ini'

        # Create a new directory with .mat file in it to compute the hybrid parameters
        cwd = os.getcwd()
        hybrid_dir = 'C:\\Temp\\hybrid'
        i = 0
        while os.path.isdir(hybrid_dir):
            i += 1
            hybrid_dir = 'C:\\Temp\\hybrid_{0:1d}'.format(i)

        test_cases = hybrid_dir + '\\testcases'
        os.makedirs(test_cases + '\\CT')

        self.set_option('Phase','orientation','file','Fiber')
        self.set_option('Phase','orientation_file','specimen.dof','Fiber')
        self.set_option('Phase','orientation_format','Digimat','Fiber')

        self.set_option('Analysis','failure_stop_analysis','off')

        self.print_daf( test_cases + '\\CT\\' + self.daf_name )
        shutil.copy(dirpath + '\\specimen.dof',test_cases + '\\CT')


        params = hybrid_dir + '\\hybrid_param'
        os.makedirs(params)

        subp = Popen([Hybrid_Command, exec_name, config_file, test_cases, '7', params], stdout=PIPE, stderr=PIPE)
        out, err = subp.communicate()

        self.read_daf(params + '\\CT\\' + self.mat_name, savefile=False)

        self.delete_option('Phase','orientation','Fiber')
        self.delete_option('Phase','orientation_file','Fiber')
        self.delete_option('Phase','orientation_format','Fiber')

        self.set_option('Analysis','failure_stop_analysis','off')
        self.set_option('Analysis','stiffness_reduction','on')
        self.set_option('Analysis','stiffness_reduction_threshold_deletion','0.99')
        self.set_option('Analysis','stiffness_reduction_max_eq_strain','0.2')

        self.print_daf(self.mat_name)
        self.print_daf(self.daf_name)
        shutil.rmtree(hybrid_dir)


############################################################
#  Object for a Woven Material
############################################################
class Woven(Daf):

    ############################################################
    #  Run the digimat program to create a new fabx file for the analysis
    ############################################################
    def create_fabric(self):
        """ Run the Digimat .daf file to create the .fabx file in MF-GUI then kill GUI """
        fabFile = self.get_option('woven','fabric_file').strip().replace('"','')

        pwd = os.getcwd()
        dafAbs = os.path.abspath(self.daf_name)
        subp = Popen([FABcommand,'input='+dafAbs], stdout=PIPE, stderr=PIPE)

        while True:
            next_line = subp.stdout.readline()
            if next_line == '' and subp.poll() is not None:
                break
            sys.stdout.flush()
            if 'fabx' in next_line:
                tmp_fab_file = next_line.replace('written: ','').strip()
                sys.stdout.flush()
                os.system('taskkill /im digimatGUI.exe /f')

        if tmp_fab_file == self.fabric_name:
            pass
        else:
            shutil.copy(tmp_fab_file,self.fabric_name)
        self.update_fabric_daf()
        return tmp_fab_file

    def update_fabric_daf(self):
        """ Update the fiber volume fractions in the phases of the
        microstructure for the computed FVF from the .fabx file """
        fab_tree = ET.parse(self.fabric_name)
        root = fab_tree.getroot()
        # Compute the total volume in the microstructure
        x = float(root.find("./Fabric/Data/X").text)
        y = float(root.find("./Fabric/Data/Y").text)
        z = float(root.find("./Fabric/Data/Z").text)
        total_volume = x*y*z

        # Iterate over the yarns in the fabric
        yarns = root.findall("./Fabric/Yarns/YarnInFabric/*")
        yarn_lengths = root.findall("./Fabric/Parameters/YarnLength/*")
        lengths = { l.tag : float(l.text) for l in yarn_lengths }

        warp = 0.
        weft = 0.
        for yarn in yarns:
            d1 = float(yarn.find("AveD1").text)
            d2 = float(yarn.find("AveD2").text)
            l = lengths[yarn.tag]
            vol =  pi*d1*d2*l/4.

            if 'warp' in yarn.find("StructurePosition").text:
                warp += vol
                self.warp_area = vol/l
            else:
                weft += vol
                self.weft_area = vol/l

        self.warp_vf = warp/total_volume
        self.weft_vf = weft/total_volume
        self.matrix_vf = 1.0 - self.warp_vf - self.weft_vf
        print('Warp:', self.warp_vf)
        print('Weft:', self.weft_vf)
        print('Matrix:', self.matrix_vf)

        self.set_option('phase','volume_fraction',
                        '{:12.6e}'.format(self.warp_vf),'WarpYarn')
        self.set_option('phase','volume_fraction',
                        '{:12.6e}'.format(self.weft_vf),'WeftYarn')
        self.set_option('phase','volume_fraction',
                        '{:12.6e}'.format(self.matrix_vf),'Matrix')

    ############################################################
    #  Run through the commands to make and analyze an analysis with fabric
    ############################################################
    def processWoven(self):
        """ Process the given daf file by creating the fabric file then running the
        analysis and reading the eng results """
        # 1. Read the daf file
        print('# 1. Read the daf file')
        self.read_daf()

        # 2. Run the daf file to create the fabric
        print('# 2. Run the daf file to create the fabric')
        fab_file_name = self.create_fabric()
        #cwd = os.getcwd()
        #shutil.copy(fabFileName,'.')

        # 3. Change the edaf file to use the temporary fabric file
        print('# 3. Change the .daf file to use the temporary fabric file')
        #this_daf.set_option('woven','fabric_file',fabFileName)
        self.print_daf()
        shutil.copy(self.daf_name,self.mat_name)

        # 4. Run the MF analysis
        print('# 4. Run the MF analysis')
        self.run_MF()

        # 5. Read the eng file
        print('# 5. Read the eng file')
        results = self.read_eng()
        return results

############################################################
#  Object for an SFRP material
############################################################
class SFRP(Daf):
    """ This class contains the functions for manipulating a discontinuous
    fiber material """

    def __init__(self,job_name='tmp',mat_file = None):
        """ Creates a material file off the baseline SFRP .daf file to modify based
        on the given material, access the methods in the .daf file with mySFRP.daf"""
        if mat_file is None:
            mat_file = 'SFRP_Baseline.mat'
            shutil.copy(dirpath+'/'+mat_file, mat_file)
        #super(SFRP, self).__init__(job_name, mat_file) # Call __init__ from the base class
        super().__init__(job_name, mat_file) # Call __init__ from the base class
        self.rename(job_name)

    def add_failure_envelope(self, npoints = 20):
        """ Add the call to compute failure envelope to the mat file """
        # Check the ANALYSIS section for 'DPP' keywords
        FE_options = { 'DPP_method': 'failure_envelope',
                       'DPP_nb_divisions': str(int(npoints)),
                       'DPP_temperature': '0.000000000000000e+000',
                       'DPP_load_type': '11_vs_22',
                       'DPP_strain_based_outputs': 'off',
                       'DPP_outputs': 'Custom,failure' }
        sec_num = self.find_section('ANALYSIS')

        for key, value in FE_options.items():
            self.options[sec_num].append(key)
            self.values[sec_num].append(value)

    ############################################################
    #  Analyze the performance of the material
    ############################################################
    def strength_vs_angle(self, n=10):
        """ Compute the strength of the material vs loading angle at n points """
        angles = np.linspace(0., 90., n)
        strengths = []
        strains = []
        for i in range(n):
            self.set_load_orientation(angles[i])
            self.run_analysis()
            self.get_stress_strain()
            strengths.append(self.maxs)
            strains.append(self.maxe)
        self.strengths = strengths
        self.strains = strains
        self.angles = angles
        return strengths, strains, angles

    def failure_envelope(self, n=10, shear = 0.):
        """ Compute the failure of the material at various loading combinations
        shear is a scalar from 0 to 1 indicating the presence of shear vs normal loads"""
        angles = np.linspace(0, 2.*pi, n+1)
        s11s = []
        s22s = []
        s12s = []
        self.run_analysis()
        self.read_eng()
        for i in range(n):
            # Compute the stresses in the one and two direction for the given direction
            s11 = cos(angles[i])*(1-shear)
            s22 = sin(angles[i])*(1-shear)
            s12 = shear
            e = np.dot(self.results['C'],np.array([s11,s22,0.,s12,0.,0.]))
            e = e*(0.1/np.linalg.norm(e))
            self.loading_2D(e[0],e[1],e[3])
            self.run_analysis()
            self.get_stress_strain()
            stress = self.s[-1]
            s11s.append(stress[0])
            s22s.append(stress[1])
            s12s.append(stress[3])
        self.s11s = s11s
        self.s22s = s22s
        self.s12s = s12s

        return s11s, s22s, s12s

    def calibrate_matrix_E(self,E11,v12, mass_fraction=None, density=1.6, orientation=[0.5, 0.5, 0]):
        """ optimize the matrix modulus to match the nominal composite stiffness
        of the SFRP system (2D Random) """
        if not with_scipy:
            print('scipy was not loaded, calibration not possible')
            return None
        self.set_fiber(orientation=orientation)
        if mass_fraction == None:
            mass_fraction = self.get_fiber_mass_fraction()
        else:
            self.set_fiber(mass_frac=mass_fraction)

        if density == None:
            density = self.get_composite_density()
        else:
            self.set_matrix_density(density)

        # Define the starting point for calibration
        Em0 = self.get_matrix('modulus')
        v0 = self.get_matrix('poisson')

        # optimization function
        def func(x):
            self.set_matrix(modulus=x[0]*Em0, poisson=x[1]*v0)
            results = self.run_analysis()
            print('Matrix Modulus = {0:>8.1f}MPa , Composite Modulus = {1:>8.1f}MPa'.format(x[0] * Em0, results['Ex']))
            print('Matrix poisson = {0:>8.4f}    , Composite Poisson = {1:>8.4f}'.format(x[1] * v0, results['vxy']))
            return ((results['Ex']-E11)/E11)**2 + ((results['vxy']-v12)/v12)**2

        # Run the optimization and return the optimized values
        results = minimize(func, [1., 1.], method='BFGS')
        x = results['x']
        Em = x[0]*Em0
        vm = x[1]*v0
        #self.set_matrix(modulus=Em, poisson=vm)
        return results

    def compute_composite_strength(self):
        """ Analyzes the current material model for the effective composite strengths

        """
        # Axial Tension
        self.set_axial_loading(e=0.05)

############################################################
#  Utility Functions
############################################################
def _strip_data(engFileLines):
    """ Parse the lines from the eng file to get the numerical data """
    data = []
    for line in engFileLines.splitlines():
        try:
            lineData = [ float(x) for x in line.split()]
            data.append(lineData)
        except:
            pass
    return np.array(data)
