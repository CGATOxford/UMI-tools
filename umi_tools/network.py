'''
network.py - Network methods for dealing with UMIs
=========================================================

:Author: Tom Smith
:Release: $Id$
:Date: |today|
:Tags: Python UMI

'''
import collections
import itertools
import sys
import numpy as np

import pyximport
pyximport.install(build_in_temp=False)

try:
    from umi_tools._dedup_umi import edit_distance
    import umi_tools.Utilities as U

except:
    from _dedup_umi import edit_distance
    import Utilities as U

sys.setrecursionlimit(10000)


def breadth_first_search(node, adj_list):
    searched = set()
    found = set()
    queue = set()
    queue.update((node,))
    found.update((node,))

    while len(queue) > 0:
        node = (list(queue))[0]
        found.update(adj_list[node])
        queue.update(adj_list[node])
        searched.update((node,))
        queue.difference_update(searched)

    return found


def recursive_search(node, adj_list):
    children = adj_list[node]
    children = [x for x in children if x not in recursive_search.component]
    for child in children:
        recursive_search.component.update((child,))
        recursive_search.component.update(
            recursive_search(child, adj_list))
    return recursive_search.component


def breadth_first_search_recursive(node, adj_list):
    try:
        recursive_search.component = set((node,))
        return recursive_search(node, adj_list)

    except RecursionError as error:
        U.info('Recursion Error: %s' % error)
        return breadth_first_search(node, adj_list)


def remove_umis(adj_list, cluster, nodes):
    '''removes the specified nodes from the cluster and returns
    the remaining nodes '''

    # list incomprehension: for x in nodes: for node in adj_list[x]: yield node
    nodes_to_remove = set([node
                           for x in nodes
                           for node in adj_list[x]] + nodes)

    return cluster - nodes_to_remove


class ReadClusterer:
    '''A functor that clusters a bundle of reads.

    Optionally:

      - identify the parent UMIs and return:
         - selected reads
         - umis
         - counts

    The initiation of the functor defines the methods:

      ** get_adj_list ** - returns the edges connecting the UMIs

      ** connected_components ** - returns clusters of connected components
                                   using the edges in the adjacency list

      ** get_best ** - returns the parent UMI(s) in the connected_components

      ** reduce_clusters ** - loops through the connected components in a
                              cluster and returns the unique reads. Optionally
                              returns lists of umis and counts per umi also

    Note: The get_adj_list and connected_components methods are not required by
    all custering methods. Where there are not required, the methods return
    None or the input parameters.

    '''

    # "get_best" methods #

    def _get_best_min_account(self, cluster, adj_list, counts):
        ''' return the min UMI(s) need to account for cluster'''
        if len(cluster) == 1:
            return list(cluster)

        sorted_nodes = sorted(cluster, key=lambda x: counts[x],
                              reverse=True)

        for i in range(len(sorted_nodes) - 1):
            if len(remove_umis(adj_list, cluster, sorted_nodes[:i+1])) == 0:
                return sorted_nodes[:i+1]

    def _get_best_higher_counts(self, cluster, counts):
        ''' return the UMI with the highest counts'''
        if len(cluster) == 1:
            return list(cluster)[0]
        else:
            sorted_nodes = sorted(cluster, key=lambda x: counts[x],
                                  reverse=True)
            return sorted_nodes[0]

    def _get_best_percentile(self, cluster, counts):
        ''' return all UMIs with counts >1% of the
        median counts in the cluster '''

        if len(cluster) == 1:
            return list(cluster)
        else:
            threshold = np.median(list(counts.values()))/100
            return [read for read in cluster if counts[read] > threshold]

    def _get_best_null(self, cluster, counts):
        ''' return all UMIs in the cluster'''

        return list(cluster)

    # "get_adj_list" methods #

    def _get_adj_list_adjacency(self, umis, counts, threshold):
        ''' identify all umis within hamming distance threshold'''

        adj_list = {umi: [] for umi in umis}
        for umi1, umi2 in itertools.combinations(umis, 2):
            if edit_distance(umi1, umi2) <= threshold:
                adj_list[umi1].append(umi2)
                adj_list[umi2].append(umi1)

        return adj_list

    def _get_adj_list_directional(self, umis, counts, threshold=1):
        ''' identify all umis within the hamming distance threshold
        and where the counts of the first umi is > (2 * second umi counts)-1'''

        adj_list = {umi: [] for umi in umis}
        for umi1, umi2 in itertools.combinations(umis, 2):
            if edit_distance(umi1, umi2) <= threshold:
                if counts[umi1] >= (counts[umi2]*2)-1:
                    adj_list[umi1].append(umi2)
                if counts[umi2] >= (counts[umi1]*2)-1:
                    adj_list[umi2].append(umi1)

        return adj_list

    def _get_adj_list_null(self, umis, counts, threshold):
        ''' for methods which don't use a adjacency dictionary'''
        return None

    # "get_connected_components" methods #

    def _get_connected_components_adjacency(self, umis, graph, counts):
        ''' find the connected UMIs within an adjacency dictionary'''

        # TS: TO DO: Work out why recursive function does lead to same
        # final output. Then uncomment below

        #if len(graph) < 10000:
        #    self.search = breadth_first_search_recursive
        #else:
        #    self.search = breadth_first_search

        found = set()
        components = list()

        for node in sorted(graph, key=lambda x: counts[x], reverse=True):
            if node not in found:
                #component = self.search(node, graph)
                component = breadth_first_search(node, graph)
                found.update(component)
                components.append(component)

        return components

    def _get_connected_components_null(self, umis, adj_list, counts):
        ''' for methods which don't use a adjacency dictionary'''
        return umis

    # "group" methods #

    def _group_unique(self, clusters, adj_list, counts):
        ''' return groups for unique method'''
        if len(clusters) == 1:
            groups = [clusters]
        else:
            groups = [[x] for x in clusters]

        return groups

    def _group_directional(self, clusters, adj_list, counts):
        ''' return groups for directional method'''

        observed = set()
        groups = []

        for cluster in clusters:
            if len(cluster) == 1:
                groups.append(list(cluster))
                observed.update(cluster)
            else:
                # need to remove any node which has already been observed
                temp_cluster = []
                for node in cluster:
                    if node not in observed:
                        temp_cluster.append(node)
                        observed.add(node)
                groups.append(temp_cluster)

        return groups

    def _group_adjacency(self, clusters, adj_list, counts):
        ''' return groups for adjacency method'''

        groups = []

        for cluster in clusters:
            if len(cluster) == 1:
                groups.append(list(cluster))

            else:
                observed = set()
                sorted_nodes = sorted(
                    cluster, key=lambda x: counts[x], reverse=True)
                temp_groups = []

                for i in range(len(sorted_nodes) - 1):
                    top_node = sorted_nodes[i]
                    observed.add(top_node)
                    latest_temp_group = [top_node]

                    for connected_node in adj_list[top_node]:
                        if connected_node == top_node:
                            continue
                        elif connected_node in observed:
                            continue
                        else:
                            latest_temp_group.append(connected_node)
                            observed.add(connected_node)

                    # need to remove the top node from any
                    # other group where it is found
                    new_temp_groups = []
                    for temp_group in temp_groups:
                        temp_group = [x for x in temp_group if x != top_node]
                        new_temp_groups.append(temp_group)
                    temp_groups = new_temp_groups

                    temp_groups.append(latest_temp_group)

                    if len(remove_umis(adj_list, cluster, sorted_nodes[:i+1])) == 0:
                        groups.extend(temp_groups)
                        break

        return groups

    def _group_cluster(self, clusters, adj_list, counts):
        ''' return groups for cluster or directional methods'''

        groups = []
        for cluster in clusters:
            groups.append(sorted(cluster, key=lambda x: counts[x], reverse=True))

        return groups

    # "reduce_clusters" methods #

    def _reduce_clusters_multiple(self, bundle, clusters,
                                  adj_list, counts, stats=False):
        ''' collapse clusters down to the UMI(s) which account for the cluster
        using the adjacency dictionary and return the list of final UMIs'''

        # TS - the "adjacency" variant of this function requires an adjacency
        # list to identify the best umi, whereas the other variants don't
        # As temporary solution, pass adj_list to all variants
        reads = []
        final_umis = []
        umi_counts = []

        for cluster in clusters:
            parent_umis = self.get_best(cluster, adj_list, counts)
            reads.extend([bundle[umi]["read"] for umi in parent_umis])

            if stats:
                final_umis.extend(parent_umis)
                umi_counts.extend([counts[umi] for umi in parent_umis])

        return reads, final_umis, umi_counts

    def _reduce_clusters_single(self, bundle, clusters,
                                adj_list, counts, stats=False):
        ''' collapse clusters down to the UMI which accounts for the cluster
        using the adjacency dictionary and return the list of final UMIs'''

        reads = []
        final_umis = []
        umi_counts = []

        for cluster in clusters:
            parent_umi = self.get_best(cluster, counts)
            reads.append(bundle[parent_umi]["read"])

            if stats:
                final_umis.append(parent_umi)
                umi_counts.append(sum([counts[x] for x in cluster]))

        return reads, final_umis, umi_counts

    def _reduce_clusters_no_network(self, bundle, clusters,
                                    adj_list, counts, stats=False):
        ''' collapse down to the UMIs which accounts for the cluster
        and return the list of final UMIs'''

        reads = []
        final_umis = []
        umi_counts = []
        parent_umis = self.get_best(clusters, counts)
        reads.extend([bundle[umi]["read"] for umi in parent_umis])

        if stats:
            final_umis.extend(parent_umis)
            umi_counts.extend([counts[umi] for umi in parent_umis])

        return reads, final_umis, umi_counts

    def __init__(self, cluster_method="directional"):
        ''' select the required class methods for the cluster_method'''

        if cluster_method == "adjacency":
            self.get_adj_list = self._get_adj_list_adjacency
            self.get_connected_components = self._get_connected_components_adjacency
            self.get_best = self._get_best_min_account
            self.reduce_clusters = self._reduce_clusters_multiple
            self.get_groups = self._group_adjacency

        elif cluster_method == "directional":
            self.get_adj_list = self._get_adj_list_directional
            self.get_connected_components = self._get_connected_components_adjacency
            self.get_best = self._get_best_higher_counts
            self.reduce_clusters = self._reduce_clusters_single
            self.get_groups = self._group_directional

        elif cluster_method == "cluster":
            self.get_adj_list = self._get_adj_list_adjacency
            self.get_connected_components = self._get_connected_components_adjacency
            self.get_best = self._get_best_higher_counts
            self.reduce_clusters = self._reduce_clusters_single
            self.get_groups = self._group_cluster

        elif cluster_method == "percentile":
            self.get_adj_list = self._get_adj_list_null
            self.get_connected_components = self._get_connected_components_null
            self.get_best = self._get_best_percentile
            self.reduce_clusters = self._reduce_clusters_no_network
            # percentile method incompatible with defining UMI groups
            self.get_groups = None

        if cluster_method == "unique":
            self.get_adj_list = self._get_adj_list_null
            self.get_connected_components = self._get_connected_components_null
            self.get_best = self._get_best_null
            self.reduce_clusters = self._reduce_clusters_no_network
            self.get_groups = self._group_unique

    def __call__(self, bundle, threshold, stats=False, further_stats=False,
                 deduplicate=True):
        ''' '''

        umis = bundle.keys()

        len_umis = [len(x) for x in umis]
        assert max(len_umis) == min(len_umis), (
            "not all umis are the same length(!):  %d - %d" % (
                min(len_umis), max(len_umis)))

        counts = {umi: bundle[umi]["count"] for umi in umis}

        adj_list = self.get_adj_list(umis, counts, threshold)

        clusters = self.get_connected_components(umis, adj_list, counts)

        if not deduplicate:
            groups = [list(x) for x in
                      self.get_groups(clusters, adj_list, counts)]
            return bundle, groups, counts

        reads, final_umis, umi_counts = self.reduce_clusters(
            bundle, clusters, adj_list, counts, stats)

        if further_stats:
            topologies = collections.Counter()
            nodes = collections.Counter()

            if len(clusters) == len(umis):
                topologies["single node"] = len(umis)
                nodes[1] = len(umis)
            else:
                for cluster in clusters:
                    if len(cluster) == 1:
                        topologies["single node"] += 1
                        nodes[1] += 1
                    else:
                        most_con = max([len(adj_list[umi]) for umi in cluster])

                        if most_con == len(cluster):
                            topologies["single hub"] += 1
                            nodes[len(cluster)] += 1
                        else:
                            topologies["complex"] += 1
                            nodes[len(cluster)] += 1

        else:
            topologies = None
            nodes = None

        return reads, final_umis, umi_counts, topologies, nodes
