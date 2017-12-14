# Sophisticated Single Cell Model Builder
# By: Brenyn Jungmann
# Version: 1.0
# Notes: This is a program to create a Sophisticated Single Cell Model without having to write any code.  There are
# multiple steps to do this, see more details below.
# 1. Program to read in a HOC Cell Template from NEURON's Cell Builder and convert it to a Python Cell Template for use
#    with NetPyNe and Python.
# 2. Program to take in user input from a Python configuration file and create a runnable model without the user
#    having to write any code.
# import ConfigParser


# ----- Start of Step 1: Template Converter -----
# Function to check if a variable is a number.
# Code used from https://www.pythoncentral.io/how-to-check-if-a-string-is-a-number-in-python-including-unicode/
def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        pass

    try:
        import unicodedata
        unicodedata.numeric(s)
        return True
    except (TypeError, ValueError):
        pass

    return False


# This function determines how I connect different sections of the cell
def find_connect_vars(start, current_line):
    for a in range(start, len(current_line)):
        stripped_line = current_line[a].strip(',')
        current_line[a] = stripped_line
        first_pareth = current_line[a].find("(")
        sec_pareth = current_line[a].find(")")
        if a == 1 or a == 6:
            if first_pareth + 1 == sec_pareth - 1:
                first_connect_loc = current_line[a][first_pareth + 1]
            else:
                first_connect_loc = current_line[a][first_pareth + 1:sec_pareth]
        elif a == 2 or a == 7:
            if first_pareth + 1 == sec_pareth - 1:
                sec_connect_loc = current_line[a][first_pareth + 1]
            else:
                sec_connect_loc = current_line[a][first_pareth + 1:sec_pareth]
        current_line[a] = current_line[a][0:first_pareth]
    return current_line, first_connect_loc, sec_connect_loc


# This function creates a string that for the basic shape of the cell.
def basic_shape_func(cur_section, sp_line, sec_names, sec_names_flag, complexity_flag):
    section_confirmed = 0
    count_var = 0
    basic_shape_string = ''
    if complexity_flag == 0:
        for b in range(0, len(sp_line)):
            var = ""
            if sp_line[b] == sec_names[sec_names_flag]:
                section_confirmed = 1
            elif section_confirmed:
                split = sp_line[b].strip('{},()')
                sp_line[b] = split
                split = sp_line[b].strip('{,()')
                sp_line[b] = split
                while len(sp_line[b]) > 0 and not sp_line[b][-1].isalpha():
                    if sp_line[b][-1] != "(":
                        var = sp_line[b][-1] + var
                        sp_line[b] = sp_line[b][:-1]
                if sp_line[b] == "pt3dclear":
                    basic_shape_string += "{2}{2}h.{0}(sec=self.{1})\n".format(sp_line[b], sec_names[sec_names_flag],
                                                                               tab)
                elif sp_line[b] == "pt3dadd":
                    basic_shape_string += "{2}{2}h.{0}({1}, ".format(sp_line[b], var, tab)
                else:
                    if count_var == 2:
                        basic_shape_string += "{0}, sec=self.{1})\n".format(var, sec_names[sec_names_flag])
                        count_var = 0
                    else:
                        basic_shape_string += "{0}, ".format(var)
                        count_var += 1
    else:
        if len(sp_line) == 2:
            cur_section = sp_line[0]
            if cur_section in section_names:
                if cur_section[-1] != "]" and is_sec_names_array[section_names.index(cur_section)]:
                    cur_section = cur_section + '[0]'
            for b in range(0, len(sp_line)):
                split = sp_line[b].strip('{},()')
                sp_line[b] = split
                split = sp_line[b].strip('{,()')
                sp_line[b] = split
                if sp_line[b] == "pt3dclear":
                    basic_shape_string += "{2}{2}h.{0}(sec=self.{1})\n".format(sp_line[b], cur_section, tab)
        elif len(sp_line) == 4:
            for b in range(0, len(sp_line)):
                var = ""
                split = sp_line[b].strip('{},()')
                sp_line[b] = split
                split = sp_line[b].strip('{,()')
                sp_line[b] = split
                while len(sp_line[b]) > 0 and not sp_line[b][-1].isalpha():
                    var = sp_line[b][-1] + var
                    sp_line[b] = sp_line[b][:-1]
                if sp_line[b] == "pt3dadd":
                    basic_shape_string += "{2}{2}h.{0}{1}, ".format(sp_line[b], var, tab)
                elif sp_line[b] == "pt3dstyle":
                    basic_shape_string += "{2}{2}h.{0}{1}, ".format(sp_line[b], var, tab)
                else:
                    if count_var == 2:
                        basic_shape_string += "{0}, sec=self.{1})\n".format(var, cur_section)
                        count_var = 0
                    else:
                        basic_shape_string += "{0}, ".format(var)
                        count_var += 1
    return basic_shape_string, cur_section


# Define list variables
python_file_string = "from neuron import h\n\n\n# Define a Class for Cell\nclass Cell(object):\n    # Create an " \
                     "initialization function that runs when a cell is created\n    # Create sections of the cell, " \
                     "connect cell sections, and create functions to initialize each section\n    def __init__(self):\n"
section_names = []
sec_list_names = []
functionList = []
sec_param_string = []
sec_param_num = []
is_sec_names_array = []

# Define string variables
init_sec_function = ""
shape_string = ""
basic_shape_func_list = []
subset_string_array = []
geo_string = ""
geom_nseg_string = ""
biophysics_string = ""
tab = "    "
variable = ""
current_function = ""
current_section = ''
bio_cur_section = ''

# Define flag variables
topology_flag = -1
basic_shape_flag = -1
bshape_list_flag = -1
complex_bshape_flag = 0
complex_basic_flag = -1
subset_flag = -1
geo_flag = -1
geo_flag_complex = -1
geom_nseg_flag = -1
biophysics_flag = -1
bio_flag_complex = -1
bio_forsec_flag = -1
section_names_flag = -1
insert_flag = 0
count = 0
sec_confirmed = 0
parentheses_count = 0
is_sec_array_count = 0

# Open HOC Template
hoc_temp = open('Cell_Template.hoc', 'r')

# Create and Open New Python Template File
python_temp = open('Cell_Template.py', 'w')

# Read HOC template line by line to determine if the file is complex or not.
# A complex template is determined by what is called in the basic_shape() function in the hoc template.
for line in hoc_temp:
    split_line = line.split()
    if split_line:
        if split_line[0] == "}" and bshape_list_flag > -1:
            bshape_list_flag = -1
        elif bshape_list_flag > -1:
            if len(split_line) == 1:
                basic_shape_func_list.append(str(split_line))
                complex_bshape_flag = 1
        if len(split_line) > 2:
            # If I find a line with basic_shape() in it I check for the data inside the procedure and then
            # create a list accordingly.
            if split_line[1] == "basic_shape()":
                bshape_list_flag = 0

# Start Reading file from beginning.
hoc_temp.seek(0, 0)

# Read HOC template line by line looking for information that we need to convert to Python
for line in hoc_temp:
    split_line = line.split()
    if split_line:
        # If I find a line starting with create I split, strip, and append the line in proper format to the python
        # template string
        if split_line[0] == "}":
            if topology_flag > -1:
                topology_flag = -1
            if subset_flag > -1:
                subset_flag = -1
            if geom_nseg_flag > -1:
                geom_nseg_flag = -1
            if geo_flag_complex > -1:
                geo_flag_complex = -1
        if complex_basic_flag > -1 and parentheses_count == 0:
            complex_basic_flag = -1
        if split_line[0] == "create":
            for i in range(1, len(split_line)):
                s = split_line[i].strip(',')
                split_line[i] = s
                if split_line[i][len(split_line[i]) - 1] == ']':
                    length = len(split_line[i]) - 1
                    first_bracket = split_line[i].find("[")
                    size_of_array = split_line[i][first_bracket + 1:length]
                    split_line[i] = split_line[i][0:first_bracket]
                    section_names.append(split_line[i])
                    is_sec_names_array.append(1)
                    create_string = "{2}{2}self.{0} = [h.Section(name='{0}[%d]' % i, cell=self) for i in " \
                                    "xrange({1})]\n".format(split_line[i], size_of_array, tab)
                else:
                    section_names.append(split_line[i])
                    is_sec_names_array.append(0)
                    is_sec_array_count += 1
                    create_string = "{1}{1}self.{0} = h.Section(name='{0}', cell=self)\n".format(split_line[i], tab)
                python_file_string += create_string
        # If I find a line starting with connect I split it and append the line in proper format to the python
        # template string
        elif topology_flag > -1:
            topology_flag += 1
            if split_line[0] == "connect":
                split_line, first_conn_loc, sec_conn_loc = find_connect_vars(1, split_line)
                connect_var1 = split_line[1][0:]
                connect_var2 = split_line[2]
                for i in range(0, len(section_names)):
                    if is_sec_names_array[i] == 1:
                        if connect_var1 == section_names[i]:
                            connect_var1 = connect_var1 + '[0]'
                        if connect_var2 == section_names[i]:
                            connect_var2 = connect_var2 + '[0]'
                connect_string = "{4}{4}self.{0}.connect(self.{1}({2}), {3})\n" \
                    .format(connect_var1, connect_var2, sec_conn_loc, first_conn_loc, tab)
                python_file_string += connect_string
            elif split_line[0] == "for":
                loop_lower = split_line[3].strip(',')
                loop_upper = split_line[4]
                split_line, first_conn_loc, sec_conn_loc = find_connect_vars(6, split_line)
                connect_string = "{0}{0}for i in range({1}, {2}):\n{0}{0}{0}self.{4}.connect(self.{3}({5}), {6})\n" \
                    .format(tab, loop_lower, int(loop_upper) + 1, split_line[6].strip(','), split_line[7],
                            first_conn_loc,
                            sec_conn_loc)
                python_file_string += connect_string
        elif -1 < basic_shape_flag < len(section_names):
            if complex_bshape_flag == 0:
                section_names_flag += 1
                shape_string = basic_shape_func(current_section, line, section_names, section_names_flag, complex_bshape_flag)
                basic_shape_flag += 1
                sec_param_string.append(shape_string)
                sec_param_num.append(section_names_flag)
        elif complex_basic_flag > -1:
            if len(split_line) > 1:
                if '{' in split_line[1]:
                    parentheses_count += 1
            elif split_line[0] == '}':
                parentheses_count -= 1
            shape_string, current_section = basic_shape_func(current_section, split_line, section_names,
                                                             section_names_flag, complex_bshape_flag)
            sec_param_string.append(shape_string)
        elif -1 < subset_flag:
            if subset_flag == 0:
                for i in range(0, len(split_line)):
                    s = split_line[i].strip(',')
                    split_line[i] = s
                    if split_line[i] != "objref":
                        sec_list_names.append(split_line[i])
                        subset_string_array.append("{0}{0}self.{1} = h.SectionList()\n".format(tab, split_line[i]))
            elif subset_flag > 0 and len(split_line) == 2:
                for j in range(0, len(section_names)):
                    if split_line[0] == section_names[j]:
                        for k in range(0, len(sec_list_names)):
                            if split_line[1] == "{0}.append()".format(sec_list_names[k]):
                                s = split_line[1].strip(')')
                                split_line[1] = s
                                subset_string_array.append("{0}{0}self.{1}sec=self.{2})\n".format(tab, split_line[1],
                                                                                                  split_line[0]))
            elif subset_flag > 0 and len(split_line) > 2:
                if split_line[0] == "for":
                    loop_lower = split_line[1].strip(',')
                    loop_lower = loop_lower[-1:]
                    s = split_line[4].strip(')')
                    split_line[4] = s
                    subset_string_array.append(
                        "{0}{0}for i in range({1}, {2}):\n{0}{0}{0}self.{4}sec=self.{3})\n".format(
                            tab, loop_lower, split_line[2], split_line[3], split_line[4]))
            subset_flag += 1
        # elif -1 < geo_flag < len(section_names):
        #     # print "Geometry: " + str(section_names_flag)
        #     section_names_flag += 1
        #     for i in range(1, len(split_line)):
        #         if not (split_line[i] == "{" or split_line[i] == "}" or split_line[i] == "="):
        #             if is_number(split_line[i]):
        #                 geo_string += "= {0}\n".format(split_line[i])
        #                 sec_param_string.append(geo_string)
        #                 sec_param_num.append(section_names_flag)
        #             else:
        #                 geo_string = "{2}{2}{0}.{1} ".format(section_names[geo_flag], split_line[i], tab)
        #     geo_flag += 1
        elif -1 < geo_flag_complex:
            if split_line[0] == 'forsec':
                geo_string += '{0}{0}for sec in self.{1}:\n{0}{0}{0}sec.{2} = {3}\n{0}{0}{0}sec.{4} = {5}\n'\
                    .format(tab, split_line[1], split_line[3], split_line[5], split_line[6], split_line[8])
            else:
                geo_string += '{0}{0}self.{1}.{2} = {3}\n{0}{0}self.{1}.{4} = {5}\n'\
                    .format(tab, split_line[0], split_line[2], split_line[4], split_line[5], split_line[7])
        elif -1 < geom_nseg_flag < len(section_names):
            if complex_bshape_flag == 0:
                section_names_flag += 1
                for i in range(1, len(split_line)):
                    if not (split_line[i] == "{" or split_line[i] == "}" or split_line[i] == "="):
                        if is_number(split_line[i]):
                            geom_nseg_string += "= {0}\n".format(split_line[i])
                            sec_param_string.append(geom_nseg_string)
                            sec_param_num.append(section_names_flag)
                        else:
                            geom_nseg_string = "{2}{2}{0}.{1} ".format(section_names[geom_nseg_flag], split_line[i], tab)
                geom_nseg_flag += 1
        # Simple hoc template biophysics data converter
        elif -1 < biophysics_flag < len(section_names):
            for i in range(0, len(split_line)):
                if split_line[i] == "}":
                    biophysics_flag += 1
                    sec_confirmed = 0
                elif sec_confirmed:
                    if not split_line[i] == "=":
                        if is_number(split_line[i]):
                            biophysics_string += " = {0}\n".format(split_line[i])
                            sec_param_string.append(biophysics_string)
                            sec_param_num.append(section_names_flag)
                        elif insert_flag == 1:
                            biophysics_string = "{2}{2}{0}.insert('{1}')\n".format(section_names[biophysics_flag],
                                                                                   split_line[i], tab)
                            sec_param_string.append(biophysics_string)
                            sec_param_num.append(section_names_flag)
                            insert_flag = 0
                        else:
                            if split_line[i] == "insert":
                                insert_flag = 1
                            else:
                                biophysics_string = "{2}{2}{0}.{1}".format(section_names[biophysics_flag], split_line[i]
                                                                           , tab)
                elif split_line[i] == section_names[biophysics_flag]:
                    sec_confirmed = 1
                    section_names_flag += 1
        # If the file is complex (Came from NeuroMorpho.com) then I process the data differently. This is for all the
        # biophysics data in the hoc file.
        elif -1 < bio_flag_complex:
            if bio_forsec_flag < 0:
                if bio_forsec_flag == -1 and split_line[0] == '}':
                    bio_flag_complex = -1
                elif split_line[0] == 'forsec':
                    bio_cur_section = split_line[1]
                    bio_forsec_flag = 1
                    sec_confirmed = 1
                    biophysics_string += "{1}{1}for sec in self.{0}:\n".format(bio_cur_section, tab)
                    # bio_forsec_flag = 1
                else:
                    if split_line[0] in section_names:
                        bio_cur_section = split_line[0]
                        bio_forsec_flag = 0
                        sec_confirmed = 1
            elif bio_forsec_flag == 0:
                for i in range(0, len(split_line)):
                    if split_line[i] == "}":
                        bio_forsec_flag = -1
                    else:
                        if not split_line[i] == "=":
                            if is_number(split_line[i]):
                                biophysics_string += " = {0}\n".format(split_line[i])
                            elif insert_flag == 1:
                                biophysics_string += "{2}{2}self.{0}.insert('{1}')\n".format(bio_cur_section,
                                                                                                split_line[i], tab)
                                insert_flag = 0
                            else:
                                if split_line[i] == "insert":
                                    insert_flag = 1
                                else:
                                    biophysics_string += "{2}{2}self.{0}.{1}".format(bio_cur_section, split_line[i],
                                                                                        tab)
            elif bio_forsec_flag == 1:
                for i in range(0, len(split_line)):
                    if split_line[i] == "}":
                        bio_forsec_flag = -1
                    else:
                        if not split_line[i] == "=":
                            if is_number(split_line[i]):
                                biophysics_string += " = {0}\n".format(split_line[i])
                            elif insert_flag == 1:
                                biophysics_string += "{2}{2}{2}sec.insert('{1}')\n".format(bio_cur_section,
                                                                                           split_line[i], tab)
                                insert_flag = 0
                            else:
                                if split_line[i] == "insert":
                                    insert_flag = 1
                                else:
                                    biophysics_string += "{2}{2}{2}sec.{1}".format(bio_cur_section, split_line[i],
                                                                                   tab)
        # If I find "proc" then I set the appropriate flag that corresponds to the hoc function name and then next
        # time through the program will move into those conditional statements to read and format the data properly
        elif split_line[0] == "proc":
            if split_line[1] in str(basic_shape_func_list):
                section_names_flag = 0
                complex_basic_flag = 0
                parentheses_count = 1
                current_function = str(split_line[1])
            elif split_line[1] == "topol()":
                topology_flag = 0
            elif split_line[1] == "basic_shape()":
                if complex_bshape_flag == 0:
                    section_names_flag = -1
                    basic_shape_flag = 0
                    count = 0
            elif split_line[1] == "subsets()":
                section_names_flag = -1
                subset_flag = 0
            elif split_line[1] == "geom_nseg()":
                section_names_flag = -1
                geom_nseg_flag = 0
            elif split_line[1] == "geom()":
                if complex_bshape_flag == 0:
                    section_names_flag = -1
                    geo_flag = 0
                elif complex_bshape_flag == 1:
                    section_names_flag = -1
                    geo_flag_complex = 0
            elif split_line[1] == "biophys()":
                if complex_bshape_flag == 0:
                    section_names_flag = -1
                    biophysics_flag = 0
                    sec_confirmed = 0
                elif complex_bshape_flag == 1:
                    section_names_flag = -1
                    bio_flag_complex = 0

# Once you get to this point everything is read in properly from the file and saved in the designated variables
# that are needed to be written to the output template file. From here I start writing all that information to a
# variable python_file_string that will be wrote one time to a file and will create a python template
# Depending on the complexity of the hoc template being converted a different format may be used in the python
# template file but that does not effect the Cell that will be created from it.

# Initialize the section lists of the Cell
for i in range(0, len(sec_list_names)):
    python_file_string += "{0}".format(subset_string_array[i])
python_file_string += "{0}{0}self.init_section_list()\n".format(tab)
python_file_string += "{0}{0}self.init_cell_basic_shape()\n".format(tab)
python_file_string += "{0}{0}self.init_geometry()\n".format(tab)
python_file_string += "{0}{0}self.init_biophysics()\n\n".format(tab)

if complex_bshape_flag == 0:
    for i in range(0, len(section_names)):
        func_name = "{1}# Function that defines the section {0}'s geometry and biophysics\n{1}" \
                    "def init_{0}(self):\n{1}{1}{0} = self.{0}\n\n".format(section_names[i], tab)
        functionList.append(func_name)
        init_sec_function = "{1}{1}self.init_{0}()\n".format(section_names[i], tab)
        python_file_string += init_sec_function
    python_file_string += "\n"


python_file_string += "{0}# Function to define the Cell's basic shape\n{0}def init_cell_basic_shape(self):\n" \
    .format(tab)

if complex_bshape_flag == 0:
    for i in range(0, len(section_names)):
        for j in range(0, len(sec_param_num)):
            if sec_param_num[j] == i:
                python_file_string += sec_param_string[j]
        python_file_string += "\n"
else:
    for j in range(0, len(sec_param_string)):
        python_file_string += sec_param_string[j]

python_file_string += "\n{0}# Function to define all section lists\n{0}def init_section_list(self):\n".format(tab)
for i in range(len(sec_list_names), len(subset_string_array)):
    python_file_string += "{0}".format(subset_string_array[i])
python_file_string += "\n"

python_file_string += '{0}# Function to define all sections geometry\n{0}def init_geometry(self):\n'\
                       .format(tab) + geo_string + '\n'

if complex_bshape_flag == 0:
    for i in range(0, len(section_names)):
        python_file_string += functionList[i]
else:
    python_file_string += '{0}# Function to define all sections biophysics\n{0}def init_biophysics(self):\n' \
                          .format(tab) + biophysics_string

# Write to file and close hoc and python template files.
python_temp.write(python_file_string)
hoc_temp.close()
python_temp.close()

# ----- End of Step 1: Template Converter -----

quit()
