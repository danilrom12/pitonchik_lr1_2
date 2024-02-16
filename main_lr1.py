import numpy as np
import re

def parse_file():
  # Local variables
  fl_str_in_trans = 0
  fl_col_out_cap = 0
  fl_gate_start = 0
  fl_gate_cap = 0
  fl_gate_fall = 0
  fl_gate_rise = 0
  fl_gate_fall_trans = 0
  fl_gate_rise_trans = 0

  # Table headers
  in_trans_list = []
  out_cap_list = []

  # Gate analysis
  fall_tab = []
  rise_tab = []
  fall_trans_tab = []
  rise_trans_tab = []

  # Parse file
  with open('lab1.lib', 'r', encoding='utf-8') as infile:
    for line in infile:
      # Fetch transition time
      if (fl_str_in_trans):
        in_trans_list = np.array(line.replace("\n", "").split(","), dtype='float')
        fl_str_in_trans = 0

      # Fetch output capacitance
      if (fl_col_out_cap):
        out_cap_list = np.array(line.replace("\n", "").split(","), dtype='float')
        fl_col_out_cap = 0

      # Fetch gate capacitance
      if (fl_gate_cap):
        gate_cap = float(line)
        fl_gate_cap = 0

      # Fetch gate fall
      if (fl_gate_fall == 1):
        if (line != '\n'):
          fall_tab.append(np.array(line.replace("\n", "").split(","), dtype='float'))
        else:
          fl_gate_fall += 1

      # Fetch gate rise
      if (fl_gate_rise == 1):
        if (line != '\n'):
          rise_tab.append(np.array(line.replace("\n", "").split(","), dtype='float'))
        else:
          fl_gate_rise += 1

      # Fetch gate fall transition
      if (fl_gate_fall_trans == 1):
        if (line != '\n'):
          fall_trans_tab.append(np.array(line.replace("\n", "").split(","), dtype='float'))
        else:
          fl_gate_fall_trans += 1

      # Fetch gate rise transition
      if (fl_gate_rise_trans == 1):
        if (line != '\n'):
          rise_trans_tab.append(np.array(line.replace("\n", "").split(","), dtype='float'))
        else:
          fl_gate_rise_trans += 1

      # Check gate EOL
      if (fl_gate_start and fl_gate_fall == 2 and fl_gate_rise == 2 and fl_gate_fall_trans == 2 and fl_gate_rise_trans == 2):
        fl_gate_start = 0
        fl_gate_cap = 0
        fl_gate_fall = 0
        fl_gate_rise = 0
        fl_gate_fall_trans = 0
        fl_gate_rise_trans = 0

      # Check input transition
      if (re.match('string_input_transtion', line)):
        fl_str_in_trans = 1
      # Check output capacitance
      if (re.match('column_output_capacitance', line)):
        fl_col_out_cap = 1
      # Check gate tables start
      if (re.match(gate_name, line)):
        fl_gate_start = 1
      # Check gate capacitance
      if (fl_gate_start and re.match('capacitance', line)):
        fl_gate_cap = 1
      # Check gate fall
      if (fl_gate_start and re.match('cell_fall', line)):
        fl_gate_fall = 1
      # Check gate rise
      if (fl_gate_start and re.match('cell_rise', line)):
        fl_gate_rise = 1
      # Check gate fall transition
      if (fl_gate_start and re.match('fall_transition', line)):
        fl_gate_fall_trans = 1
      # Check gate rise transition
      if (fl_gate_start and re.match('rise_transition', line)):
        fl_gate_rise_trans = 1
  return in_trans_list, out_cap_list, fall_tab, rise_tab, fall_trans_tab, rise_trans_tab

def find_borders(in_data, in_data_list, in_closest_index):
  # Get left and right index
  if (in_data >= in_data_list[in_closest_index]):
    ind_l = in_closest_index
    ind_r = in_closest_index + 1
  else:
    ind_l = in_closest_index - 1
    ind_r = in_closest_index
  return ind_l, ind_r

def find_borders_value(i_l, i_r, list):
  return list[i_l], list[i_r]

def bilinear_interpol(table, str_l, str_r, col_l, col_r, hdr_str_l, hdr_str_r, hdr_col_l, hdr_col_r, in_str, in_col):
  aprox_str_col_l = ((hdr_col_r - in_col) / (hdr_col_r - hdr_col_l) * table[str_l][col_l] +
               (in_col - hdr_col_l) / (hdr_col_r - hdr_col_l) * table[str_l][col_r])
  aprox_str_col_r = ((hdr_col_r - in_col) / (hdr_col_r - hdr_col_l) * table[str_r][col_l] +
               (in_col - hdr_col_l) / (hdr_col_r - hdr_col_l) * table[str_r][col_r])
  aprox_result = ((hdr_str_r - in_str) / (hdr_str_r - hdr_str_l) * aprox_str_col_l +
                  (in_str - hdr_str_l) / (hdr_str_r - hdr_str_l) * aprox_str_col_r)
  return aprox_result

# Input data
gate_name = "BUF"
gate_trans_type = "n"
gate_in_trans = "17"
gate_out_cap = "23"
# Convert input
gate_in_trans = float(gate_in_trans)
gate_out_cap = float(gate_out_cap)

# Parse file
in_trans_list, out_cap_list, fall_tab, rise_tab, fall_trans_tab, rise_trans_tab = parse_file()
# Check out of range input
if ((gate_in_trans < in_trans_list[0] or gate_in_trans > in_trans_list[-1]) or
  (gate_out_cap < out_cap_list[0] or gate_out_cap > out_cap_list[-1]) or
  not((gate_trans_type == "p") or (gate_trans_type == "n"))):
  print("Out of range")
else:
  # Get closest table index
  in_trans_ind = np.abs(in_trans_list - gate_in_trans).argmin()
  out_cap_ind = np.abs(out_cap_list - gate_out_cap).argmin()
  # Find table border
  tr_i_l, tr_i_r = find_borders(gate_in_trans, in_trans_list, in_trans_ind)
  cap_i_l, cap_i_r = find_borders(gate_out_cap, out_cap_list, out_cap_ind)
  tr_v_l, tr_v_r = find_borders_value(tr_i_l, tr_i_r, in_trans_list)
  cap_v_l, cap_v_r = find_borders_value(cap_i_l, cap_i_r, out_cap_list)
  # Interpolate value
  if (gate_trans_type == "p"):
    print("Closest value")
    print("Rise:", rise_tab[in_trans_ind][out_cap_ind], "ps")
    print("Rise transition:", rise_trans_tab[in_trans_ind][out_cap_ind], "ps")
    print("Bilinear interpolation")
    print("Rise:", bilinear_interpol(rise_tab, tr_i_l, tr_i_r, cap_i_l, cap_i_r, tr_v_l, tr_v_r, cap_v_l, cap_v_r, gate_in_trans, gate_out_cap), "ps")
    print("Rise:", bilinear_interpol(rise_trans_tab, tr_i_l, tr_i_r, cap_i_l, cap_i_r, tr_v_l, tr_v_r, cap_v_l, cap_v_r, gate_in_trans, gate_out_cap), "ps")
  else:
    print("Closest value")
    print("Fall:", fall_tab[in_trans_ind][out_cap_ind], "ps")
    print("Fall transition:", fall_trans_tab[in_trans_ind][out_cap_ind], "ps")
    print("Bilinear interpolation")
    print("Fall:", bilinear_interpol(fall_tab, tr_i_l, tr_i_r, cap_i_l, cap_i_r, tr_v_l, tr_v_r, cap_v_l, cap_v_r, gate_in_trans, gate_out_cap), "ps")
    print("Fall transition:", bilinear_interpol(fall_trans_tab, tr_i_l, tr_i_r, cap_i_l, cap_i_r, tr_v_l, tr_v_r, cap_v_l, cap_v_r, gate_in_trans, gate_out_cap), "ps")