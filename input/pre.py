# def main():
#     # input_fname = 'Cluster_output_normalized_IL_Acc95_2_75'
#     input_fname = 'Subcluster_output_for_cluster_400'
#     output_fname = input_fname + '.in'
#     fp = open(input_fname)
#     lines = fp.readlines()
#     fp.close()

#     clusters = {}
#     for i, line in enumerate(lines):
#         if i == 0:
#             continue
#         # motif_id, cluster_no, fam_id = line.strip().split('\t')
#         # cluster_id = 'cluster_' + str(cluster_no)


#         motif_id, cluster_no, subcluster_no, u_subcluster_id, subcluster_no, fam_id = line.strip().split('\t')
#         cluster_id = 'subcluster_' + str(subcluster_no)

#         pieces = motif_id.strip().split('_')
#         loop = pieces[1] + '_' + pieces[2] + ':' + '_'.join(pieces[3:])
#         if cluster_id not in clusters:
#             clusters[cluster_id] = []
#         clusters[cluster_id].append(loop)

#     print(clusters)
#     fp = open(output_fname, 'w')
#     for cluster_id in sorted(clusters):
#         # print(cluster_id, end='')
#         fp.write(cluster_id)
#         for loop in clusters[cluster_id]:
#             # print(',' + loop, end='')
#             fp.write(',' + loop)
#         # print('')
#         fp.write('\n')
#     fp.close()

def main():
    # input_fname = 'Cluster_output_normalized_IL_Acc95_2_75'
    # input_fname = 'Subcluster_output_for_cluster_400'
    input_fname = 'cluster_output.csv'
    output_fname = input_fname + '.in'
    fp = open(input_fname)
    lines = fp.readlines()
    fp.close()

    clusters = {}
    for i, line in enumerate(lines):
        if i == 0:
            continue
        # motif_id, cluster_no, fam_id = line.strip().split('\t')
        # cluster_id = 'cluster_' + str(cluster_no)


        # motif_id, cluster_no, subcluster_no, u_subcluster_id, subcluster_no, fam_id = line.strip().split('\t')
        motif_id, cluster_no, subcluster_no, fam_id = line.strip().split('\t')
        cluster_id = 'subcluster_' + str(subcluster_no)

        # pieces = motif_id.strip().split('_')
        # pieces = motif_id.replace(':', '_').strip().split('_')
        # loop = pieces[1] + '_' + pieces[2] + ':' + '_'.join(pieces[3:])
        loop = motif_id
        if cluster_id not in clusters:
            clusters[cluster_id] = []
        clusters[cluster_id].append(loop)

    print(clusters)
    fp = open(output_fname, 'w')
    for cluster_id in sorted(clusters):
        # print(cluster_id, end='')
        fp.write(cluster_id)
        for loop in clusters[cluster_id]:
            # print(',' + loop, end='')
            loop = '_'.join(loop.strip().split('_')[:2]) + ':' + '_'.join(loop.strip().split('_')[2:])
            fp.write(',' + loop)
        # print('')
        fp.write('\n')
    fp.close()
            




if __name__ == '__main__':
    main()