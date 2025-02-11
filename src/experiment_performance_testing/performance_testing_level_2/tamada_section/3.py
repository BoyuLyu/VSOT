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
import multiprocessing

def make_plane(mesh, path, vertices):
    """
    Compute cross-sections along a path on a 3D mesh.
    """
    print(f"make_plane called with {len(vertices)} vertices")
    out_slice_3D = []
    out_area = []
    for ind, p0 in enumerate(vertices):
        # Compute normal vectors and cross-section
        # Your existing logic here...
        pass
    return out_area, out_slice_3D

def process_file(all_files2, rootFolder, outputFolder, namex):
    try:
        print(f"Processing {namex}")
        namex = str(namex.split('.')[0])
        mesh = trimesh.load_mesh(os.path.join(rootFolder, namex+ '.off'), process=False)
        # path_vertices = mesh.vertices[:10]  # Replace with actual path logic
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
        pathx = nx.shortest_path(g,
                                source=start,
                                target=end,
                                weight='length')
        path_vertices = mesh.vertices[sspath]
        # generate a series of cross sections along the path
        area_list, slice_3d = make_plane(mesh, pathx, path_vertices)
        area_list, slice_3d = make_plane(mesh, [], path_vertices)
        print(f"Completed make_plane for {namex}")
    except Exception as e:
        print(f"Error in process_file: {e}")

if __name__ == '__main__':
    rootFolder = "/work/boyu/EM_astrocyte/materials_for_paper_VSOT/VSOT_data_paired_with_codes/data/performance_testing_level_2_segmentation/large_dataset_3_annotator_L2345_w_stubby/gt_final_800/stubby_samples/off_file"
    outputFolder = "/work/boyu/EM_astrocyte/materials_for_paper_VSOT/VSOT_data_paired_with_codes/data/performance_testing_level_2_segmentation/large_dataset_3_annotator_L2345_w_stubby/segmentation_results/tamada_result_stubby"
    all_files = os.listdir(rootFolder)
    all_files2 = os.listdir(outputFolder)

    with multiprocessing.Pool(2) as p:
        p.starmap(process_file, [(all_files2, rootFolder, outputFolder, namex) for namex in all_files])

