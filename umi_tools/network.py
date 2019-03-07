'''
network.py - Network methods for dealing with UMIs
=========================================================

'''

from __future__ import absolute_import
import collections
import itertools
import sys
import regex
import numpy as np

from umi_tools._dedup_umi import edit_distance
import umi_tools.Utilities as U
import umi_tools.whitelist_methods as whitelist_methods

sys.setrecursionlimit(10000)


def breadth_first_search(node, adj_list):
    searched = set()
    queue = set()
    queue.update((node,))
    searched.update((node,))

    while len(queue) > 0:
        node = queue.pop()
        for next_node in adj_list[node]:
            if next_node not in searched:
                queue.update((next_node,))
                searched.update((next_node,))

    return searched


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


def get_substr_slices(umi_length, idx_size):
    '''
    Create slices to split a UMI into approximately equal size substrings
    Returns a list of tuples that can be passed to slice function
    '''
    cs, r = divmod(umi_length, idx_size)
    sub_sizes = [cs + 1] * r + [cs] * (idx_size - r)
    offset = 0
    slices = []
    for s in sub_sizes:
        slices.append((offset, offset + s))
        offset += s
    return slices


def build_substr_idx(umis, umi_length, min_edit):
    '''
    Build a dictionary of nearest neighbours using substrings, can be used
    to reduce the number of pairwise comparisons.
    '''
    substr_idx = collections.defaultdict(
        lambda: collections.defaultdict(set))
    slices = get_substr_slices(umi_length, min_edit + 1)
    for idx in slices:
        for u in umis:
            u_sub = u[slice(*idx)]
            substr_idx[idx][u_sub].add(u)
    return substr_idx


def iter_nearest_neighbours(umis, substr_idx):
    '''
    Added by Matt 06/05/17
    use substring dict to get (approximately) all the nearest neighbours to
    each in a set of umis.
    '''
    for i, u in enumerate(umis, 1):
        neighbours = set()
        for idx, substr_map in substr_idx.items():
            u_sub = u[slice(*idx)]
            neighbours = neighbours.union(substr_map[u_sub])
        neighbours.difference_update(umis[:i])
        for nbr in neighbours:
            yield u, nbr


class UMIClusterer:
    '''A functor that clusters a dictionary of UMIs and their counts.
    The primary return value is either a list of representative UMIs
    or a list of lists where each inner list represents the contents of
    one cluster.

    Optionally:

      - identify the parent UMIs and return:
         - selected reads
         - umis
         - counts

    The initiation of the functor defines the methods:

      ** get_adj_list ** - returns the edges connecting the UMIs

      ** get_connected_components ** - returns clusters of connected components
                                       using the edges in the adjacency list

      ** get_groups ** - returns the groups of umis,
                         with the parent umi at position 0

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

    def _get_best_percentile(self, cluster, counts):
        ''' return all UMIs with counts >1% of the
        median counts in the cluster '''

        if len(cluster) == 1:
            return list(cluster)
        else:
            threshold = np.median(list(counts.values()))/100
            return [read for read in cluster if counts[read] > threshold]

    # "get_adj_list" methods #

    def _get_adj_list_adjacency(self, umis, counts, threshold):
        ''' identify all umis within hamming distance threshold'''

        adj_list = {umi: [] for umi in umis}
        if len(umis) > 25:
            umi_length = len(umis[0])
            substr_idx = build_substr_idx(umis, umi_length, threshold)
            iter_umi_pairs = iter_nearest_neighbours(umis, substr_idx)
        else:
            iter_umi_pairs = itertools.combinations(umis, 2)
        for umi1, umi2 in iter_umi_pairs:
            if edit_distance(umi1, umi2) <= threshold:
                adj_list[umi1].append(umi2)
                adj_list[umi2].append(umi1)

        return adj_list

    def _get_adj_list_directional(self, umis, counts, threshold=1):
        ''' identify all umis within the hamming distance threshold
        and where the counts of the first umi is > (2 * second umi counts)-1'''

        adj_list = {umi: [] for umi in umis}
        if len(umis) > 25:
            umi_length = len(umis[0])
            substr_idx = build_substr_idx(umis, umi_length, threshold)
            iter_umi_pairs = iter_nearest_neighbours(umis, substr_idx)
        else:
            iter_umi_pairs = itertools.combinations(umis, 2)
        for umi1, umi2 in iter_umi_pairs:
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

        # TS: TO DO: Work out why recursive function doesn't lead to same
        # final output. Then uncomment below

        # if len(graph) < 10000:
        #    self.search = breadth_first_search_recursive
        # else:
        #    self.search = breadth_first_search

        found = set()
        components = list()

        for node in sorted(graph, key=lambda x: counts[x], reverse=True):
            if node not in found:
                # component = self.search(node, graph)
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
                cluster = sorted(cluster, key=lambda x: counts[x],
                                 reverse=True)
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

                lead_umis = self._get_best_min_account(cluster,
                                                       adj_list, counts)
                observed.update(lead_umis)

                for lead_umi in lead_umis:
                    connected_nodes = set(adj_list[lead_umi])
                    groups.append([lead_umi] +
                                  list(connected_nodes - observed))
                    observed.update(connected_nodes)

        return groups

    def _group_cluster(self, clusters, adj_list, counts):
        ''' return groups for cluster or directional methods'''

        groups = []
        for cluster in clusters:
            groups.append(sorted(cluster, key=lambda x: counts[x],
                                 reverse=True))

        return groups

    def _group_percentile(self, clusters, adj_list, counts):
        ''' Return "groups" for the the percentile method. Note
        that grouping isn't really compatible with the percentile
        method. This just returns the retained UMIs in a structure similar
        to other methods '''

        retained_umis = self._get_best_percentile(clusters, counts)
        groups = [[x] for x in retained_umis]

        return groups

    def __init__(self, cluster_method="directional"):
        ''' select the required class methods for the cluster_method'''

        self.max_umis_per_position = 0
        self.total_umis_per_position = 0
        self.positions = 0

        if cluster_method == "adjacency":
            self.get_adj_list = self._get_adj_list_adjacency
            self.get_connected_components = self._get_connected_components_adjacency
            self.get_groups = self._group_adjacency

        elif cluster_method == "directional":
            self.get_adj_list = self._get_adj_list_directional
            self.get_connected_components = self._get_connected_components_adjacency
            self.get_groups = self._group_directional

        elif cluster_method == "cluster":
            self.get_adj_list = self._get_adj_list_adjacency
            self.get_connected_components = self._get_connected_components_adjacency
            self.get_groups = self._group_cluster

        elif cluster_method == "percentile":
            self.get_adj_list = self._get_adj_list_null
            self.get_connected_components = self._get_connected_components_null
            # percentile method incompatible with defining UMI groups
            self.get_groups = self._group_percentile

        elif cluster_method == "unique":
            self.get_adj_list = self._get_adj_list_null
            self.get_connected_components = self._get_connected_components_null
            self.get_groups = self._group_unique

    def __call__(self, umis, threshold):
        '''umis is a dictionary that maps UMIs to their counts'''

        counts = umis
        umis = list(umis.keys())

        self.positions += 1

        number_of_umis = len(umis)

        self.total_umis_per_position += number_of_umis

        if number_of_umis > self.max_umis_per_position:
            self.max_umis_per_position = number_of_umis

        len_umis = [len(x) for x in umis]

        assert max(len_umis) == min(len_umis), (
            "not all umis are the same length(!):  %d - %d" % (
                min(len_umis), max(len_umis)))

        adj_list = self.get_adj_list(umis, counts, threshold)
        clusters = self.get_connected_components(umis, adj_list, counts)
        final_umis = [list(x) for x in
                      self.get_groups(clusters, adj_list, counts)]

        return final_umis


class ReadDeduplicator:
    '''This is a wrapper for applying the UMI methods to bundles of BAM reads.
    It is currently a pretty transparent wrapper on UMIClusterer. Basically
    taking a read bundle, extracting the UMIs and Counts, running UMIClusterer
    and returning the results along with annotated reads

    In addition, if UMI whitelist options are provided, the read
    clusters are filtered against the UMI whitelist'''

    def __init__(self, options):

        self.UMIClusterer = UMIClusterer(cluster_method=options.method)

        if options.filter_umi:
            self.umi_whitelist = whitelist_methods.getUserDefinedBarcodes(
                options.umi_whitelist,
                options.umi_whitelist_paired,
                deriveErrorCorrection=False)[0]
            self.umi_whitelist_counts = collections.Counter()

            U.info("Length of UMI whitelist: %i" % len(self.umi_whitelist))

        else:
            self.umi_whitelist = None

    def __call__(self, bundle, threshold):
        '''Process the the bundled reads according to the method specified
        in the constructor. Return signature is:

        reads, final_umis, umi_counts, topologies, nodes

        reads:        predicted best reads for deduplicated position
        final_umis:   list of predicted parent UMIs
        umi_counts:   Sum of read counts for reads represented by the
                      corresponding UMI
        '''

        umis = bundle.keys()
        counts = {umi: bundle[umi]["count"] for umi in umis}

        clusters = self.UMIClusterer(counts, threshold)

        if self.umi_whitelist:
            # check the "top" UMI is in the whitelist, if not, discard
            # the whole group

            final_umis = []
            umi_counts = []

            for cluster in clusters:
                cluster_count = sum(counts[x] for x in cluster)
                if cluster[0].decode() in self.umi_whitelist:
                    final_umis.append(cluster[0])
                    umi_counts.append(cluster_count)
                    self.umi_whitelist_counts[cluster[0].decode()] += cluster_count

                else:
                    self.umi_whitelist_counts["Non-whitelist UMI"] += cluster_count
        else:
            final_umis = [cluster[0] for cluster in clusters]
            umi_counts = [sum(counts[umi] for umi in cluster)
                          for cluster in clusters]

        reads = [bundle[umi]["read"] for umi in final_umis]

        return (reads, final_umis, umi_counts)


class CellClusterer:
    '''A functor that clusters a dictionary of cell barcodes and their counts.
    The primary return value is either a list of representative UMIs
    or a list of lists where each inner list represents the contents of
    one cluster.

    '''

    def _get_best_min_account(self, cluster, adj_list, counts):
        ''' return the min UMI(s) need to account for cluster'''
        if len(cluster) == 1:
            return list(cluster)

        sorted_nodes = sorted(cluster, key=lambda x: counts[x],
                              reverse=True)

        for i in range(len(sorted_nodes) - 1):
            if len(remove_umis(adj_list, cluster, sorted_nodes[:i+1])) == 0:
                return sorted_nodes[:i+1]

    def _get_adj_list_directional(self, umis, counts):
        ''' identify all umis within the hamming distance threshold
        and where the counts of the first umi is > (2 * second umi counts)-1'''

        adj_list = {umi: [] for umi in umis}

        if self.fuzzy_match:
            for umi1 in umis:
                # we need a second regex for some insertions,
                # e.g UMI1 = "ATCG", UMI2 = "ATTC"
                comp_regex_err = regex.compile("(%s){e<=1}" % str(umi1))
                comp_regex_del = regex.compile("(%s){i<=1}" % str(umi1)[::-1])
                for umi2 in umis:
                    if umi1 == umi2:
                        continue
                    if counts[umi1] >= (counts[umi2]*self.dir_threshold):
                        if (max(len(umi1), len(umi2)) -
                            min(len(umi1), len(umi2))) > 1:
                            continue
                        if (comp_regex_err.match(str(umi2)) or
                            comp_regex_del.match(str(umi2))):
                            adj_list[umi1].append(umi2)
        else:
            for umi1, umi2 in itertools.combinations(umis, 2):
                if edit_distance(umi1, umi2) <= 1:
                    if counts[umi1] >= (counts[umi2]*2)-1:
                        adj_list[umi1].append(umi2)
                    if counts[umi2] >= (counts[umi1]*2)-1:
                        adj_list[umi2].append(umi1)

        return adj_list

    def _get_connected_components_adjacency(self, graph, counts):
        ''' find the connected UMIs within an adjacency dictionary'''

        found = set()
        components = list()

        for node in sorted(graph, key=lambda x: counts[x], reverse=True):
            if node not in found:
                # component = self.search(node, graph)
                component = breadth_first_search(node, graph)
                found.update(component)
                components.append(component)
        return components

    def __init__(self, cluster_method="directional",
                 dir_threshold=10, fuzzy_match=True):
        ''' select the required class methods for the cluster_method
        and set the attributes used to refine the directional method'''

        self.dir_threshold = dir_threshold
        self.fuzzy_match = fuzzy_match

        if cluster_method == "directional":
            self.get_adj_list = self._get_adj_list_directional
            self.get_connected_components = self._get_connected_components_adjacency
        else:
            raise ValueError("CellClusterer currently only supports the directional method")

    def __call__(self, umis, counts):
        '''Counts is a directionary that maps UMIs to their counts'''

        len_umis = [len(x) for x in umis]
        if not max(len_umis) == min(len_umis):
            U.warn("not all umis are the same length(!):  %d - %d" % (
                min(len_umis), max(len_umis)))

        adj_list = self.get_adj_list(umis, counts)

        clusters = self.get_connected_components(umis, adj_list, counts)

        final_umis = [list(x) for x in
                      self.get_groups(clusters, adj_list, counts)]

        return final_umis
