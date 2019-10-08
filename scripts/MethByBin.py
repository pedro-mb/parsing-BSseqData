#!/usr/bin/python


import tabix
import argparse


def get_methylation(tb_query, report, sample, stt, end, strand, out_dict, flank, bin_number, bin_number_fl, do_flanks,
                    max_depth, min_depth, feature):

    for i in tb_query:
        m_scaf, pos, m_strand, c_n, t_n, m_context, seq = i
        depth = int(t_n) + int(c_n)
        if min_depth <= depth <= max_depth:
            unit = (int(end) - int(stt) + 1) / (bin_number - 0.01)
            if strand == "+":
                if int(pos) < stt:
                    fl_dist = stt - int(pos)
                    if not do_flanks:
                        keys = 1
                        keys = report + "\t" + sample + "\t" + feature + "\t" + str(keys)
                    elif fl_dist == flank:
                        keys = -int(fl_dist / bin_number_fl) + 1
                        keys = report + "\t" + sample + "\t" + "Upstream" + "\t" + str(keys)
                    else:
                        keys = -int(fl_dist / bin_number_fl)
                        keys = report + "\t" + sample + "\t" + "Upstream" + "\t" + str(keys)
                elif stt <= int(pos) < end:
                    keys = int((int(pos) - stt + 1) / unit) + 1
                    keys = report + "\t" + sample + "\t" + feature + "\t" + str(keys)
                else:
                    fl_dist = int(pos) - end
                    if not do_flanks:
                        keys = bin_number
                        keys = report + "\t" + sample + "\t" + feature + "\t" + str(keys)
                    elif fl_dist == flank:
                        keys = int(fl_dist / bin_number_fl) + bin_number - 1
                        keys = report + "\t" + sample + "\t" + "Downstream" + "\t" + str(keys)
                    else:
                        keys = int(fl_dist/bin_number_fl) + bin_number + 1
                        keys = report + "\t" + sample + "\t" + "Downstream" + "\t" + str(keys)

            else:
                if int(pos) <= stt:
                    fl_dist = stt - int(pos)
                    if not do_flanks:
                        keys = bin_number
                        keys = report + "\t" + sample + "\t" + feature + "\t" + str(keys)
                    elif fl_dist == flank:
                        keys = int(fl_dist / bin_number_fl) + bin_number - 1
                        keys = report + "\t" + sample + "\t" + "Downstream" + "\t" + str(keys)
                    else:
                        keys = int(fl_dist / bin_number_fl) + bin_number + 1
                        keys = report + "\t" + sample + "\t" + "Downstream" + "\t" + str(keys)
                elif stt < int(pos) <= end:
                    keys = int((int(end) - int(pos) + 1) / unit) + 1
                    keys = report + "\t" + sample + "\t" + feature + "\t" + str(keys)
                else:
                    fl_dist = int(pos) - end
                    if not do_flanks:
                        keys = 1
                        keys = report + "\t" + sample + "\t" + feature + "\t" + str(keys)
                    elif fl_dist == flank:
                        keys = -int(fl_dist / bin_number_fl) + 1
                        keys = report + "\t" + sample + "\t" + "Upstream" + "\t" + str(keys)
                    else:
                        keys = -int(fl_dist / bin_number_fl)
                        keys = report + "\t" + sample + "\t" + "Upstream" + "\t" + str(keys)
            if keys in out_dict:
                out_dict[keys][0] += int(c_n)
                out_dict[keys][1] += int(t_n)
            else:
                out_dict[keys] = [int(c_n), int(t_n)]
    return out_dict


def screen_coordinates(samplefile, coord, do_flanks, flank, bin_number, bin_number_flank, max_depth, min_depth,
                       filename):
    meth_outdict = {}
    outfile = open(filename, "w")
    sample = open(samplefile, "r")
    for item in sample:
        item_list = item.strip().split()
        sample_dir, sample_id = item_list
        print sample_dir
        tbx = tabix.open(sample_dir)
        flag = 0
        gff = open(coord[0], "r")
        gen_feature = coord[1]
        for line in gff:
            flag += 1
            if flag == 1000:
                print ".",
                flag = 0
            scaffold, feat, start_c, end_c, seq_strand = line.strip().split()[0:5]

            if do_flanks:
                # do the analysis for feature a flanking regions
                stt_flank = int(start_c) - flank
                if stt_flank < 0:
                    stt_flank = 1
                else:
                    stt_flank = int(start_c) - flank + 2
                end_flank = int(end_c) + flank - 2
                query = scaffold + ":" + str(stt_flank) + "-" + str(end_flank)
                meth_freq = tbx.querys(query)
                meth_outdict = get_methylation(meth_freq, sample_dir, sample_id, int(start_c), int(end_c), seq_strand,
                                               meth_outdict, flank, bin_number, bin_number_flank, do_flanks,
                                               max_depth, min_depth, gen_feature)
            elif int(end_c) - int(start_c) >= bin_number:
                # do the analysis only for features, with size higher than binNumber
                #stt_flank = int(start_c) #+ 2
                #end_flank = int(end_c) #- 2

                query = scaffold + ":" + start_c + "-" + end_c

                meth_freq = tbx.querys(query)

                for i in meth_freq:
                    m_scaf, pos, m_strand, c_n, t_n, m_context, seq = i
                    depth = int(t_n) + int(c_n)
                    if min_depth <= depth <= max_depth:
                        unit = (int(end_c) - int(start_c) + 1) / (bin_number - 0.01)
                        if seq_strand == "+":

                            if int(start_c) <= int(pos) <= int(end_c):
                                keys = int((int(pos) - int(start_c) + 1) / unit) + 1
                                keys = sample_dir + "\t" + sample_id + "\t" + gen_feature + "\t" + str(keys)
                                if keys in meth_outdict:
                                    meth_outdict[keys][0] += int(c_n)
                                    meth_outdict[keys][1] += int(t_n)
                                else:
                                    meth_outdict[keys] = [int(c_n), int(t_n)]
                        else:
                            if int(start_c) <= int(pos) <= int(end_c):
                                keys = int((int(end_c) - int(pos) + 1) / unit) + 1
                                keys = sample_dir + "\t" + sample_id + "\t" + gen_feature + "\t" + str(keys)
                                if keys in meth_outdict:
                                    meth_outdict[keys][0] += int(c_n)
                                    meth_outdict[keys][1] += int(t_n)
                                else:
                                    meth_outdict[keys] = [int(c_n), int(t_n)]
                #meth_outdict = get_methylation(meth_freq, filename, sample_id, int(start_c), int(end_c), seq_strand,
                #                               meth_outdict, 0, bin_number, bin_number_flank, do_flanks,
                #                               max_depth, min_depth, gen_feature)
        print "\n"

    outfile.write("file\tsample\tfeature\tbin\tC_number\tT_number\tmethylation_level\n")
    for entry in sorted(meth_outdict):
        c_num, t_num = meth_outdict[entry]
        meth_freq = float(c_num) / float(c_num + t_num)
        outfile.write(entry + "\t" + str(c_num) + "\t" + str(t_num) + "\t" + str(meth_freq) + "\n")
    sample.close()
    outfile.close()


def main():
    parser = argparse.ArgumentParser(description='Filter fasta file with a list of IDs (one perline)')
    parser.add_argument('--featCoord', required=True, metavar='tab',
                        help='[REQUIRED] text file containing the coordinates of the feature to be considered '
                             'in the following format <chr>\t<feature>\t<stt>\t<end>\t<strand>\t<id>\t<product>\n '
                             '(e.g. output of gff2tab.pl)')
    parser.add_argument('--feature', required=True, metavar='str', type=str,
                        help='name of the feature being considered. e.g. transcript, exon, intron, ... ')
    parser.add_argument('--sample', required=True, metavar='str',
                        help='[REQUIRED] list of files and corresponding sample names as: \
                         <path/to/file.bgz>\t<sampleName>\n')
    parser.add_argument('--bin', metavar="INT", type=int, default=60,
                        help='the feature will be divided in this many bins. '
                             'If flanks = False, features with size < bin number '
                             'will not be considered. Default 60.')
    parser.add_argument('--out', metavar='STR', type=str, default="methylation_frequency.txt",
                        help='Output file name')
    parser.add_argument('--min_depth', metavar="INT", type=int, default=5,
                        help='minimum mC base coverage to be considered. Default: 5 ')
    parser.add_argument('--max_depth', metavar="INT", type=int, default=1000000,
                        help='maximum mC base coverage to be considered. Default: 1000000 ')
    parser.add_argument('-flanks', dest='flanking', action='store_true',
                        help= 'if selected, the analysis will include the regions flanking the feature (flanks = True)')
    parser.add_argument('--fl_length', metavar='INT', type=int, default=2000,
                        help='if flanks, size in bp of the flanking region to be considered. Default: 2000.')
    parser.add_argument('--bin_fl', metavar="INT", type=int, default=100,
                        help='if flanks, the flanking regions will be divided in this many bins. Default 100.')


    args = parser.parse_args()
    print args.sample
    screen_coordinates(args.sample, [args.featCoord, args.feature], args.flanking, args.fl_length, args.bin, args.bin_fl,
                       args.max_depth, args.min_depth, args.out)

if __name__ == "__main__":
    main()

#import tabix

#flank_length = 2000
#binNumber = 60
#binNumber_flank = 100

#minDepth = 5
#maxDepth = 1000000

#minLength = 300
#maxLength = 5000000

#Samplefile = [["../innerbark/BGI_ILU_DNA_BISU_INNERBARK_01_1_val_1_bismark_bt2_pe.deduplicated.CpG_report.txt.bgz", "innerB"]]#,
              #["../leaf/BGI_ILU_DNA_BISU_LEAF_01_1_val_1_bismark_bt2_pe.deduplicated.CpG_report.txt.bgz", "leaf"],
              #["../xylem/BGI_ILU_DNA_BISU_XYLEM_01_1_val_1_bismark_bt2_pe.deduplicated.CpG_report.txt.bgz", "xylem"],
              #["../phellem/BGI_ILU_DNA_BISU_PHELLEM_01_1_val_1_bismark_bt2_pe.deduplicated.CpG_report.txt.bgz", "phellem"]]
#Argscoord = ["../genomeSequence-transcriptGFFCoord.tab", "gene"]
#Argscoord = ["../genomeSequence-exonUniqueGFFCoord.tab", "exon"]
#Argscoord = ["../genomeSequence-intronUniqueGFFCoord.tab", "intron"]
#argscoord = ["./genome_transcriptGFF_test2.tab", "Body"]

#gff = open(argscoord[0], "r")
#gen_feature = argscoord[1]
#flanking = False
#outFile = "outfile_gene_test.tab"

#in function
#meth_outDict = {}
#outFile = open("outfile_intron.tab", "w")



