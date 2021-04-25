"""
main.py file that runs the unittests when the file is called or run using the python CLI.

Part 1 of HW 5 for CISC 471, Computational Biology.

By:
    - Hershil Devnani (20001045)
    - Rayan Shaikli (20059806)

There are two ways to run the program, outlined below. Both methods run the unittests.
Sample Usage:
  $ python -m unittest unittests.py
  $ python -m main main.py
"""

import unittest
from convolution_cyclopeptide_sequencing import *

""" Note: Tests take about ~ 3 minutes to fully complete execution"""


class TestProgrammingPartOne(unittest.TestCase):
    """Testing class for the required unittests and all Rosalind problems solved for the assignment
    """

    ###########################################################################
    #  Tests for Implement ConvolutionCyclopeptideSequencing (Rosalind ba4i)  #
    ###########################################################################
    def test_convolution_cyclopeptide_sequencing_positive_cyclic_scoring_a(self):
        # test using cyclic peptide scoring

        filename = "sample_data.txt"
        m, n, spectrum = parse_data(filename)

        convolutions, conv_dict = spectral_convolution(spectrum)
        convolutions_set = convert_convolution_set_to_tuple(conv_dict)
        expands = expand_set_with_ties(convolutions_set, m)
        leader = leaderboard_cyclopeptide_sequencing(expands, spectrum, n, scoring="cyclic")

        solution = "57-129-99-71-57-80"

        self.assertEqual(peptide_mass(leader), parent_mass(spectrum))  # make sure solution mass matches parent mass
        self.assertEqual(format_peptide_mass_seq(leader), solution)

    def test_convolution_cyclopeptide_sequencing_positive_cyclic_scoring_b(self):
        # test using cyclic peptide scoring

        filename = "convolution_cyclopeptide_sequencing_test_rosalind_ba4i_a.txt"
        m, n, spectrum = parse_data(filename)

        convolutions, conv_dict = spectral_convolution(spectrum)
        convolutions_set = convert_convolution_set_to_tuple(conv_dict)
        expands = expand_set_with_ties(convolutions_set, m)
        leader = leaderboard_cyclopeptide_sequencing(expands, spectrum, n, scoring="cyclic")

        solution = "128-129-131-131-114-128-131-113-114-87-186-147-87-114-114"

        self.assertEqual(peptide_mass(leader), parent_mass(spectrum))  # make sure solution mass matches parent mass
        self.assertEqual(format_peptide_mass_seq(leader), solution)

    def test_convolution_cyclopeptide_sequencing_positive_cyclic_scoring_c(self):
        # test using cyclic peptide scoring

        filename = "convolution_cyclopeptide_sequencing_test_rosalind_ba4i_b.txt"
        m, n, spectrum = parse_data(filename)

        convolutions, conv_dict = spectral_convolution(spectrum)
        convolutions_set = convert_convolution_set_to_tuple(conv_dict)
        expands = expand_set_with_ties(convolutions_set, m)
        leader = leaderboard_cyclopeptide_sequencing(expands, spectrum, n, scoring="cyclic")

        solution = "129-186-101-163-131-99-156-99-115-114-87-114-113-115"

        self.assertEqual(peptide_mass(leader), parent_mass(spectrum))  # make sure solution mass matches parent mass
        self.assertEqual(format_peptide_mass_seq(leader), solution)

    def test_convolution_cyclopeptide_sequencing_positive_cyclic_scoring_d(self):
        # test using cyclic peptide scoring

        filename = "convolution_cyclopeptide_sequencing_test_rosalind_ba4i_c.txt"
        m, n, spectrum = parse_data(filename)

        convolutions, conv_dict = spectral_convolution(spectrum)
        convolutions_set = convert_convolution_set_to_tuple(conv_dict)
        expands = expand_set_with_ties(convolutions_set, m)
        leader = leaderboard_cyclopeptide_sequencing(expands, spectrum, n, scoring="cyclic")

        solution = "129-156-87-163-113-163-103-71-113-147-186-131-147-97-101"

        self.assertEqual(peptide_mass(leader), parent_mass(spectrum))  # make sure solution mass matches parent mass
        self.assertEqual(format_peptide_mass_seq(leader), solution)

    def test_convolution_cyclopeptide_sequencing_negative_cyclic_scoring(self):
        # test using cyclic peptide scoring

        filename = "negative_sample.txt"
        m, n, spectrum = parse_data(filename)

        convolutions, conv_dict = spectral_convolution(spectrum)
        convolutions_set = convert_convolution_set_to_tuple(conv_dict)
        expands = expand_set_with_ties(convolutions_set, m)
        # leader = leaderboard_cyclopeptide_sequencing(expands, spectrum, n)
        leader = leaderboard_cyclopeptide_sequencing(expands, spectrum, n, scoring="cyclic")

        solution = "0"

        self.assertEqual(peptide_mass(leader), parent_mass(spectrum))  # make sure solution mass matches parent mass
        self.assertEqual(format_peptide_mass_seq(leader), solution)

    def test_convolution_cyclopeptide_sequencing_positive_linear_scoring_a(self):
        # test using linear peptide scoring

        filename = "sample_data.txt"
        m, n, spectrum = parse_data(filename)

        convolutions, conv_dict = spectral_convolution(spectrum)
        convolutions_set = convert_convolution_set_to_tuple(conv_dict)
        expands = expand_set_with_ties(convolutions_set, m)
        # leader = leaderboard_cyclopeptide_sequencing(expands, spectrum, n)
        leader = leaderboard_cyclopeptide_sequencing(expands, spectrum, n, scoring="linear")

        solution = "71-99-129-57-79-58"

        self.assertEqual(peptide_mass(leader), parent_mass(spectrum))  # make sure solution mass matches parent mass
        self.assertEqual(format_peptide_mass_seq(leader), solution)

    def test_convolution_cyclopeptide_sequencing_positive_linear_scoring_b(self):
        # test using linear peptide scoring

        filename = "convolution_cyclopeptide_sequencing_test_rosalind_ba4i_a.txt"
        m, n, spectrum = parse_data(filename)

        convolutions, conv_dict = spectral_convolution(spectrum)
        convolutions_set = convert_convolution_set_to_tuple(conv_dict)
        expands = expand_set_with_ties(convolutions_set, m)
        leader = leaderboard_cyclopeptide_sequencing(expands, spectrum, n, scoring="linear")

        solution = "87-114-113-131-128-114-131-131-129-128-114-114-87-147-72-114"

        self.assertEqual(peptide_mass(leader), parent_mass(spectrum))  # make sure solution mass matches parent mass
        self.assertEqual(format_peptide_mass_seq(leader), solution)

    def test_convolution_cyclopeptide_sequencing_positive_linear_scoring_c(self):
        # test using linear peptide scoring

        filename = "convolution_cyclopeptide_sequencing_test_rosalind_ba4i_b.txt"
        m, n, spectrum = parse_data(filename)

        convolutions, conv_dict = spectral_convolution(spectrum)
        convolutions_set = convert_convolution_set_to_tuple(conv_dict)
        expands = expand_set_with_ties(convolutions_set, m)
        leader = leaderboard_cyclopeptide_sequencing(expands, spectrum, n, scoring="linear")

        solution = "129-115-113-114-87-114-115-99-156-99-131-163-101-86-100"

        self.assertEqual(peptide_mass(leader), parent_mass(spectrum))  # make sure solution mass matches parent mass
        self.assertEqual(format_peptide_mass_seq(leader), solution)

    def test_convolution_cyclopeptide_sequencing_positive_linear_scoring_d(self):
        # test using linear peptide scoring

        filename = "convolution_cyclopeptide_sequencing_test_rosalind_ba4i_c.txt"
        m, n, spectrum = parse_data(filename)

        convolutions, conv_dict = spectral_convolution(spectrum)
        convolutions_set = convert_convolution_set_to_tuple(conv_dict)
        expands = expand_set_with_ties(convolutions_set, m)
        leader = leaderboard_cyclopeptide_sequencing(expands, spectrum, n, scoring="linear")

        solution = "103-71-113-147-186-131-147-97-101-129-156-87-78-85-113-78-85"

        self.assertEqual(peptide_mass(leader), parent_mass(spectrum))  # make sure solution mass matches parent mass
        self.assertEqual(format_peptide_mass_seq(leader), solution)

    def test_convolution_cyclopeptide_sequencing_negative_linear_scoring(self):
        # test using linear peptide scoring

        filename = "negative_sample.txt"
        m, n, spectrum = parse_data(filename)

        convolutions, conv_dict = spectral_convolution(spectrum)
        convolutions_set = convert_convolution_set_to_tuple(conv_dict)
        expands = expand_set_with_ties(convolutions_set, m)
        leader = leaderboard_cyclopeptide_sequencing(expands, spectrum, n, scoring="linear")

        solution = "0"

        self.assertEqual(peptide_mass(leader), parent_mass(spectrum))  # make sure solution mass matches parent mass
        self.assertEqual(format_peptide_mass_seq(leader), solution)

    #########################################################
    #  Tests for Convolution of a Spectrum (Rosalind ba4h)  #
    #########################################################

    def test_spectrum_convolution_a(self):
        spectrum = [0, 137, 186, 323]
        convolution, _ = spectral_convolution(spectrum)

        solution = [137, 137, 186, 186, 323, 49]

        self.assertEqual(convolution, solution)

    def test_spectrum_convolution_b(self):
        filename = "spectral_convolution_test_rosalind_ba4h.txt"
        spectrum = parse_conv_data(filename)
        convolution, _ = spectral_convolution(spectrum)

        solution = parse_conv_data("spectral_convolution_test_solution_rosalind_ba4h.txt")

        self.assertEqual(convolution, solution)

    #################################################################
    #  Tests for LeaderboardCyclopeptideSequencing (Rosalind ba4g)  #
    #################################################################

    def test_leaderboard_cyclopeptide_sequencing_a(self):
        n = 10
        spectrum = [0, 71, 113, 129, 147, 200, 218, 260, 313, 331, 347, 389, 460]
        leader_peptide = leaderboard_cyclopeptide_sequencing(DEFAULT_EXPAND_SET, spectrum, n)
        formatted_leader_peptide = format_peptide_mass_seq(leader_peptide)

        solution = [113, 147, 71, 129]
        formatted_solution = format_peptide_mass_seq(solution)

        self.assertEqual(leader_peptide, solution)
        self.assertEqual(formatted_leader_peptide, formatted_solution)

    def test_leaderboard_cyclopeptide_sequencing_b(self):
        filename = "leaderboard_sequencing_test_rosalind_ba4g.txt"
        n, spectrum = parse_leaderboard_seqeuncing_data(filename)

        leader_peptide = leaderboard_cyclopeptide_sequencing(DEFAULT_EXPAND_SET, spectrum, n, scoring="cyclic")
        formatted_leader_peptide = format_peptide_mass_seq(leader_peptide)

        solution = [115, 101, 186, 163, 163, 103, 147, 97, 115, 99, 99, 128, 137]
        formatted_solution = format_peptide_mass_seq(solution)

        self.assertEqual(leader_peptide, solution)
        self.assertEqual(formatted_leader_peptide, formatted_solution)

    ####################################################
    #  Tests for Cyclopeptide scoring (Rosalind ba4f)  #
    ####################################################

    # test scoring based on peptide string
    def test_cyclo_pepitde_scoring_a(self):
        peptide = "NQEL"
        experimental_specturm = [0, 99, 113, 114, 128, 227, 257, 299, 355, 356, 370, 371, 484]
        # peptide_score = score(peptide, experimental_specturm)
        peptide_score = score_peptide(peptide, experimental_specturm, scoring_type="cyclic")

        solution = 11

        self.assertEqual(peptide_score, solution)

        # make sure altered function based on mass seq gives same result

    def test_cyclo_pepitde_scoring_b(self):
        filename = "cyclopeptide_scoring_test_rosalind_ba4f.txt"
        peptide, experimental_specturm = parse_scoring_data(filename)
        # peptide_score = score(peptide, experimental_specturm)
        peptide_score = score_peptide(peptide, experimental_specturm, scoring_type="cyclic")

        solution = 330

        self.assertEqual(peptide_score, solution)

    # test scoring based on a given peptide mass sequence
    def test_cyclo_pepitde_mass_seq_scoring_a(self):
        peptide = [113, 129, 128, 114]
        experimental_specturm = [0, 99, 113, 114, 128, 227, 257, 299, 355, 356, 370, 371, 484]
        # peptide_score = score_mass_seq_new(peptide, experimental_specturm)
        peptide_score = score_peptide(peptide, experimental_specturm, scoring_type="cyclic")

        solution = 11

        self.assertEqual(peptide_score, solution)

    def test_cyclo_pepitde_mass_seq_scoring_b(self):
        """-----------------UPDATE----------------"""
        filename = "cyclopeptide_scoring_test_rosalind_ba4f.txt"
        peptide, experimental_specturm = parse_scoring_data(filename)
        peptide_mass_seq = peptide_to_mass_seq(peptide)
        # peptide_score = score_mass_seq_new(peptide_mass_seq, experimental_specturm)
        peptide_score = score_peptide(peptide, experimental_specturm, scoring_type="cyclic")

        solution = 330

        self.assertEqual(peptide_score, solution)

    #############################################
    #  Tests for Cyclospectrum (Rosalind ba4c)  #
    #############################################

    def test_cyclo_specturm_a(self):
        theoretical_cyclo_spectrum = cyclo_spectrum("NQEL")

        solution = [0, 113, 114, 128, 129, 227, 242, 242, 257, 355, 356, 370, 371, 484]

        self.assertEqual(theoretical_cyclo_spectrum, solution)

    def test_cyclo_specturm_b(self):
        filename = "cyclospectrum_test_rosalind_ba4c.txt"
        data = parse_spectrum_data(filename)
        theoretical_cyclo_spectrum = cyclo_spectrum(data)

        solution = [0, 57, 57, 57, 71, 71, 101, 113, 114, 128, 129, 131, 142, 147, 156, 158, 158, 186, 186, 188, 200,
                    213, 215, 215, 218, 241, 260, 269, 270, 271, 272, 289, 289, 314, 317, 326, 331, 333, 346, 371, 372,
                    383, 388, 397, 402, 403, 404, 418, 418, 427, 428, 454, 459, 475, 475, 484, 489, 500, 511, 519, 532,
                    541, 546, 549, 559, 560, 583, 590, 603, 604, 606, 612, 613, 617, 640, 647, 661, 669, 672, 674, 688,
                    697, 707, 718, 735, 759, 760, 764, 769, 789, 790, 792, 798, 800, 801, 821, 826, 830, 831, 855, 872,
                    883, 893, 902, 916, 918, 921, 929, 943, 950, 973, 977, 978, 984, 986, 987, 1000, 1007, 1030, 1031,
                    1041, 1044, 1049, 1058, 1071, 1079, 1090, 1101, 1106, 1115, 1115, 1131, 1136, 1162, 1163, 1172,
                    1172, 1186, 1187, 1188, 1193, 1202, 1207, 1218, 1219, 1244, 1257, 1259, 1264, 1273, 1276, 1301,
                    1301, 1318, 1319, 1320, 1321, 1330, 1349, 1372, 1375, 1375, 1377, 1390, 1402, 1404, 1404, 1432,
                    1432, 1434, 1443, 1448, 1459, 1461, 1462, 1476, 1477, 1489, 1519, 1519, 1533, 1533, 1533, 1590]

        self.assertEqual(theoretical_cyclo_spectrum, solution)

    def test_cyclo_spectrum_mass_seq_a(self):
        theoretical_cyclo_spectrum_mass_seq = [113, 129, 128, 114]

        theoretical_cyclo_spectrum = cyclo_spectrum_mass_seq(theoretical_cyclo_spectrum_mass_seq)

        solution = [0, 113, 114, 128, 129, 227, 242, 242, 257, 355, 356, 370, 371, 484]

        self.assertEqual(theoretical_cyclo_spectrum, solution)

    def test_cyclo_spectrum_mass_seq_b(self):
        filename = "cyclospectrum_test_rosalind_ba4c.txt"
        data = parse_spectrum_data(filename)
        theoretical_cyclo_spectrum_mass_seq = peptide_to_mass_seq(data)

        theoretical_cyclo_spectrum = cyclo_spectrum_mass_seq(theoretical_cyclo_spectrum_mass_seq)

        solution = [0, 57, 57, 57, 71, 71, 101, 113, 114, 128, 129, 131, 142, 147, 156, 158, 158, 186, 186, 188, 200,
                    213, 215, 215, 218, 241, 260, 269, 270, 271, 272, 289, 289, 314, 317, 326, 331, 333, 346, 371, 372,
                    383, 388, 397, 402, 403, 404, 418, 418, 427, 428, 454, 459, 475, 475, 484, 489, 500, 511, 519, 532,
                    541, 546, 549, 559, 560, 583, 590, 603, 604, 606, 612, 613, 617, 640, 647, 661, 669, 672, 674, 688,
                    697, 707, 718, 735, 759, 760, 764, 769, 789, 790, 792, 798, 800, 801, 821, 826, 830, 831, 855, 872,
                    883, 893, 902, 916, 918, 921, 929, 943, 950, 973, 977, 978, 984, 986, 987, 1000, 1007, 1030, 1031,
                    1041, 1044, 1049, 1058, 1071, 1079, 1090, 1101, 1106, 1115, 1115, 1131, 1136, 1162, 1163, 1172,
                    1172, 1186, 1187, 1188, 1193, 1202, 1207, 1218, 1219, 1244, 1257, 1259, 1264, 1273, 1276, 1301,
                    1301, 1318, 1319, 1320, 1321, 1330, 1349, 1372, 1375, 1375, 1377, 1390, 1402, 1404, 1404, 1432,
                    1432, 1434, 1443, 1448, 1459, 1461, 1462, 1476, 1477, 1489, 1519, 1519, 1533, 1533, 1533, 1590]

        self.assertEqual(theoretical_cyclo_spectrum, solution)

    def test_linear_specturm_a(self):
        theoretical_linear_spectrum = linear_spectrum("NQEL")

        solution = [0, 113, 114, 128, 129, 242, 242, 257, 370, 371, 484]

        self.assertEqual(theoretical_linear_spectrum, solution)

    def test_linear_specturm_b(self):
        filename = "cyclospectrum_test_rosalind_ba4c.txt"
        data = parse_spectrum_data(filename)
        theoretical_linear_spectrum = linear_spectrum(data)

        solution = [0, 57, 57, 57, 71, 71, 101, 113, 114, 128, 129, 131, 142, 147, 156, 158, 158, 186, 186, 188, 200,
                    213, 215, 215, 218, 241, 260, 269, 270, 271, 272, 289, 289, 317, 326, 331, 333, 346, 371, 372, 383,
                    388, 397, 402, 403, 404, 418, 418, 428, 454, 459, 475, 475, 484, 489, 511, 519, 532, 541, 546, 549,
                    559, 560, 590, 603, 604, 606, 612, 617, 661, 669, 672, 674, 688, 707, 735, 759, 764, 790, 792, 800,
                    801, 821, 830, 872, 893, 921, 929, 943, 950, 977, 978, 1000, 1007, 1071, 1079, 1090, 1136, 1163,
                    1193, 1218, 1276, 1349, 1404, 1462, 1590]

        self.assertEqual(theoretical_linear_spectrum, solution)

    def random(self):
        x = [0, 137, 186, 323]
        test_conv, y = spectral_convolution(x)
        print(test_conv)
        print(y)

        filename = "sample_convolution_data.txt"
        # filename = "rosalind_ba4h_8.txt"
        x = parse_conv_data(filename)
        # print(x)
        test_conv, _ = spectral_convolution(x)
        print(test_conv)
        print()
        # for i in test_conv:
        #     print(i, end=" ")
        # print()

        # filename = "sample_convolution_solution.txt"
        # solution = parse_data(filename)
        # solution.sort()
        # test_conv.sort()

        # print(test_conv == solution)

        # filename = "sample_convolution_data.txt"
        # conv_data = parse_conv_data(filename)
        # convolutions = spectral_convolution(conv_data)
        # expands = expand_set(convolutions, 3)
        # print(expands)

        # filename = "sample_scoring_data_2.txt"
        # filename = "rosalind_ba4f_2.txt"
        # peptide, spectrum = parse_scoring_data(filename)
        # print(peptide)
        # print(spectrum)
        # theoretical_spectrum = cyclo_spectrum(peptide)
        # print(theoretical_spectrum)
        # for i in theoretical_spectrum:
        #     print(i, end=" ")
        # print()
        # peptide_score = score(peptide, spectrum)
        # print(peptide_score)
        #
        # linear_spec = linear_spectrum("NQEL")
        # print(linear_spec)
        #
        # x = formatted_to_mass_seq("113-147-71-129")
        # print(x)

        # filename = "sample_scoring_data.txt"
        # peptide, spectrum = parse_scoring_data(filename)
        #
        # theo_new = linear_spectrum_mass_seq([114, 128, 129, 113])
        # print(theo_new)
        # scr = linear_score_mass_seq([114, 128, 129, 113], spectrum)

        theo_cyclo = cyclo_spectrum_mass_seq([113, 129, 128, 114])
        print(theo_cyclo)

        # filename = "leaderboard_sample.txt"
        # m, n, spectrum = parse_data(filename)
        # leader = leaderboard_cyclopeptide_sequencing(DEFAULT_EXPAND_SET, spectrum, n)
        # print(f"leader: {leader}")
        # formatted = format_peptide_mass_seq(leader)
        # print(formatted)

        # filename = "rosalind_ba4g_3.txt"
        # m, n, spectrum = parse_data(filename)
        # leader = leaderboard_cyclopeptide_sequencing(DEFAULT_EXPAND_SET, spectrum, n)
        # print(f"leader: {leader}")
        # formatted = format_peptide_mass_seq(leader)
        # print(formatted)


if __name__ == '__main__':
    unittest.main()
