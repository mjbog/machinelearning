class DigimatFile(object):
    """ Base class for all digimat files

  Initialize a DigimatFile from the file name, job name is taken by removing the 3 letter
  extension from the file.

  read_digimat_file(filename) parses the file into sections with ordered lists of
    section_names, options, and values.  Each item in the options and values lists
    are also ordered pairs.

  find_section(section,name=None) finds the list index of the given section with
    with optional name given. Search is not case-sensitive

  get_option(section,option,value=0,name=None,verbose=False) Finds the option value
    for the given section name or number and option name

  set_option(section,option_in,value=0,name=None) Sets the option in the given section
    with optional name to the given value

  delete_option(section,option,name=None) Deletes the given option name in the section
    and its corresponding value

  print_section(section, name=None) Prints all options and values for the given section
    with optional name

  print_all_sections() Prints all sections in the file with the current option, value pairs

  print_file(new_name = None, verbose=False) Prints the sections, with option value pairs
    to the given file name. If no name given, use self.file_name that was used to
    initialize the object
  """

    def __init__(self, file_name):
        self.file_name = file_name
        if file_name[-4] == ',':
            job_name = file_name[:-4]
        else:
            job_name = file_name
        self.mat_name = job_name + '.mat'
        self.daf_name = job_name + '.daf'

        self.section_names, self.options, self.values = self.read_digimat_file(file_name)

    def read_digimat_file(self, file_name):
        """ Read in the digimat file (.daf,.mat,.map) file and store the information as:
          section_names
           |-> options
           |-> values
    where section_names is an ordered list of the KEYWORD objects
    and options and values are ordered lists of lists corresponding
    to the order of the KEYWORDS in the .daf file with the options and
    values in the order they are written in the file """

        section_names = []
        all_options = []
        all_values = []

        lines = open(file_name).readlines()

        options = []
        for line in lines:
            # Skip comments
            line = line.replace('\r','')
            if line[0] is ' ' or line[0] is '#' or len(line) < 2:
                pass
            else:
                line = line[:-1]
                # Line is defining an option value pair
                if '=' in line:
                    option, value = line.split('=')
                    options.append(option.strip())
                    values.append(value.strip())
                # Line is defining a new section name
                else:
                    if not len(options) == 0:
                        # Add the previous section to the objects container
                        all_options.append(options)
                        all_values.append(values)
                    options = []
                    values = []
                    section_name = line
                    section_names.append(section_name)
                    # print 'Object Name =',sectionName
        all_options.append(options)
        all_values.append(values)

        self.section_names = section_names
        self.options = all_options
        self.values = all_values
        return section_names, all_options, all_values

    def find_section(self, section, name=None):
        """find_section(self,section,name=None)

    Find the section with the matching name if given and return its
    list index """
        i = 0
        for section_name in self.section_names:
            if not section_name.lower() == section.lower():
                i += 1
            else:
                if name is None:
                    return i
                elif not self.get_option(i, 'name').lower() == name.lower():
                    i += 1
                else:
                    return i
        return None

    def delete_section(self, section, name=None):
        """delete_section(section,name=None)

    Find the section with the given section keyword and name and delete it
    and its accompanying options and values """
        sec_no = self.find_section(section, name)

        while sec_no:
            del self.section_names[sec_no]
            del self.options[sec_no]
            del self.values[sec_no]
            sec_no = self.find_section(section, name)

    def get_option(self, section, option, value=0, name=None, verbose=False):
        """ get_option(self,section,option,value=0,name=None,verbose=False)

    Given the section name or number, return the value for the given option """
        try:
            section_number = int(section)
        except:
            section_number = self.find_section(section, name)
        try:
            options = self.options[section_number]
            values = self.values[section_number]
        except:
            if verbose:
                print('No', section, 'section named', name, 'found')
            return None

        results = []

        for i in range(len(options)):
            if options[i].lower() == option.lower():
                results.append(values[i])

        if len(results) == 1:
            return results[0]
        elif len(results) == 0:
            if verbose:
                print('No option named', option, 'found in', section, name)
            return None
        else:
            return results

    def set_option(self, section, option_in, value=0, name=None):
        """set_option(self,section,option_in,value=0,name=None)

    Replace an option in the given section with the specified value.
    If the section has multiple instances, find the section with the
    specified name """

        # print section,option,value,name
        i = self.find_section(section, name)

        j = 0
        for option in self.options[i]:
            # print option
            if not option_in.lower() == option.lower():
                j += 1
            else:
                break
        # print j,'out of',len(self.options[i])
        if j == len(self.options[i]):
            self.values[i].append(value)
            self.options[i].append(option_in)
        else:
            self.values[i][j] = value

    def delete_option(self, section, option, name=None):
        """delete_option(self,section,option,name=None)

    Delete the option in the given section with optional name """
        i = self.find_section(section, name)
        options = self.options[i]
        values = self.values[i]
        for i in reversed(range(len(options))):
            if options[i].lower() == option.lower():
                del options[i]
                del values[i]

    def print_section(self, section, name=None):
        """print_section(self,section,name=None)

    print all of the options found in that named section """
        try:
            n_sec = int(section)
            print('Section:', self.section_names[n_sec])
            options = self.options[n_sec]
            values = self.values[n_sec]
        except ValueError:
            for i in range(len(self.section_names)):
                if self.section_names[i].lower() == section.lower():
                    n_sec = i
                    opts = self.options[n_sec]
                    vals = self.values[n_sec]
                    option_dict = {opts[i]: vals[i] for i in range(len(opts))}
                    if name == None or name.lower() == option_dict['name'].lower():
                        print('Section:', self.section_names[n_sec])
                        options = opts
                        values = vals
        except IndexError:
            print('The given section,', section, ', could not be found, only',
                  len(self.section_names), 'sections are found')
            return

        try:
            for option, value in zip(options, values):
                print(option, ':', value)
        except:
            print('The given section,', section, ', could not be found')

    def print_all_sections(self):
        """ Prints all section in the digimat file using the print_section command """
        for i in range(len(self.section_names)):
            n_sec = i
            options = self.options[n_sec]
            values = self.values[n_sec]
            option_dict = {options[i]: values[i] for i in range(len(options))}
            print('Section:', self.section_names[n_sec])
            for option, value in zip(options, values):
                print('  ', option, ':', value)

    ############################################################
    #  Print the current options in the analysis file to the daf
    ############################################################
    def print_file(self, new_name=None, verbose=False, addenda=None):
        """print_file(self,new_name = None)

    Recreate the digimat file object for the analysis, if no new_name is specified,
    the function will write to the given self.file_name """
        if new_name is None:
            new_name = self.file_name

        if verbose:
            print('Write to file:', new_name)

        new_lines = ''
        for i in range(len(self.options)):

            new_lines += '\n' + 40 * '#' + '\n'
            new_lines += self.section_names[i] + '\n'
            for j in range(len(self.options[i])):
                new_lines += self.options[i][j] + ' = ' + str(self.values[i][j]) + '\n'

        if addenda:
            new_lines += addenda

        open(new_name, 'w').write(new_lines)
