# this is a reimplementation of tamada's idea based on trimesh
# original paper was based on NeuroMorph
import matplotlib
import trimesh
import numpy as np
from shapely.geometry import LineString
# import pyglet
import networkx as nx
import matplotlib.pyplot as plt
from scipy.signal import argrelmin
from scipy.signal import argrelmax
import json
import os


def make_plane(mesh, path, vertices):
    out_slice_3D = []
    out_area = []
    for ind in range(len(sspath)):
        p0 = vertices[ind]
        if ind == 0 or ind == 1:
            p_1 = vertices[1]
            p_0 = vertices[0]
            norm_m1 = p_1 - p_0
            norm_m2 = norm_m1
        else:
            pm1 = vertices[ind - 1]
            pm2 = vertices[ind - 2]
            norm_m1 = p0 - pm1
            norm_m2 = pm1 - pm2

        N = len(sspath)
        if ind == N - 1 or ind == N - 2:
            p_N = vertices[N - 1]
            p_Nm1 = vertices[N - 2]
            norm_p1 = p_N - p_Nm1
            norm_p2 = norm_p1
        else:
            pp1 = vertices[ind + 1]
            pp2 = vertices[ind + 2]
            norm_p1 = pp1 - p0
            norm_p2 = pp2 - pp1

        norm_m1 = norm_m1 / np.linalg.norm(norm_m1)
        norm_m2 = norm_m2 / np.linalg.norm(norm_m2)
        norm_p1 = norm_p1 / np.linalg.norm(norm_p1)
        norm_p2 = norm_p2 / np.linalg.norm(norm_p2)
        norm_here = (norm_m1 + norm_p1 + norm_m2 / 2 + norm_p2 / 2) / 3
        slice = mesh.section(plane_origin=p0, plane_normal=norm_here)
        if slice is not None:
            slice_2D, to_3D = slice.to_planar()
            out_slice_3D.append(slice)
            out_area.append(slice_2D.area)
    return out_area, out_slice_3D


if __name__ == '__main__':
    rootFolder = "/work/boyu/EM_astrocyte/materials_for_paper_VSOT/VSOT_data_paired_with_codes/data/performance_testing_level_2_segmentation/large_dataset_3_annotator_L2345_w_stubby/gt_final_800/uniformly_selected_samples/off_file"
    outputFolder = "/work/boyu/EM_astrocyte/materials_for_paper_VSOT/VSOT_data_paired_with_codes/data/performance_testing_level_2_segmentation/large_dataset_3_annotator_L2345_w_stubby/segmentation_results/tamada_result"
    all_files = os.listdir(rootFolder)
    all_files2 = os.listdir(outputFolder)
    for namex in all_files:
        namex = namex[:-4]
        out_name_x = namex + '.json'
        if out_name_x not in all_files2:
            print(namex)
            mesh = trimesh.load_mesh(os.path.join(rootFolder, namex + '.off'), process=False)
            # generate skeleton
            edges = mesh.edges_unique
            length = mesh.edges_unique_length
            g = nx.Graph()
            for edge, L in zip(edges, length):
                g.add_edge(*edge, length=L)
            all_pair_length = dict(nx.all_pairs_shortest_path_length(g))
            max_length = 0
            source_max_l = 0
            sink_max_l = 0
            for i in all_pair_length.keys():
                for j in all_pair_length[i].keys():
                    if all_pair_length[i][j] > max_length:
                        max_length = all_pair_length[i][j]
                        source_max_l = i
                        sink_max_l = j
            sspath = nx.shortest_path(g,
                                      source=source_max_l,
                                      target=sink_max_l,
                                      weight='length')
            start = source_max_l
            end = sink_max_l
            # run the shortest path query using length for edge weight
            path = nx.shortest_path(g,
                                    source=start,
                                    target=end,
                                    weight='length')
            path_vertices = mesh.vertices[sspath]
            # generate a series of cross sections along the path
            area_list, slice_3d = make_plane(mesh, path, path_vertices)
            area_list_array = np.array(area_list)
            if np.mean(area_list_array[0:round(len(area_list_array) / 2)]) < np.mean(
                    area_list_array[round(len(area_list_array) / 2):(len(area_list) - 1)]):
                area_list_array_diff = area_list_array[:-1] - area_list_array[1:]
                s_id = argrelmin(area_list_array_diff, order=10)
            else:
                area_list_array_diff = area_list_array[:-1] - area_list_array[1:]
                s_id = argrelmax(area_list_array_diff, order=10)
            s_id = s_id[0]
            fin_id = s_id[np.argmin(np.abs(np.array(s_id) - round(len(area_list_array) / 2)))]
            intersect_circle_vt = slice_3d[fin_id].vertices

            with open(os.path.join(outputFolder, namex + '.json'), 'w') as filehandle:
                json.dump(intersect_circle_vt.tolist(), filehandle)
