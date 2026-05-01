# MASLD META NETWORK

**MASLD META NETWORK** is an interactive R Shiny application for exploring transcriptomic meta-analysis results in metabolic dysfunction-associated steatotic liver disease (MASLD).

The application links gene-level meta-analysis results with pathway enrichment, network topology, histological disease progression, biological sex, and protein–protein interaction data. It is designed to help users move from large gene and pathway result tables toward interpretable biological modules and candidate genes associated with MASLD progression.

The app supports exploration across two histological scoring systems:

- NAFLD Activity Score / NAS
- Fibrosis stage

Using these scoring systems, MASLD META NETWORK allows users to visualize pathways, genes, protein–protein interactions, stage-dependent pathway remodeling, sex-aware enrichment patterns, and individual transcriptome-level results.

## Citation

If you use **MASLD META NETWORK** in your research, please cite the associated manuscript and this GitHub repository.

Manuscript not available yet. 

---

## Overview

MASLD progression involves coordinated changes across metabolic, inflammatory, fibrotic, and stress-response pathways. However, enrichment results are often difficult to interpret as independent tables because many biological pathways share genes and represent overlapping biological processes.

MASLD META NETWORK addresses this by representing meta-analysis-derived results as interactive networks. In these networks, nodes represent biological entities such as pathways, genes, or pathway clusters, while edges represent shared genes, cluster–gene relationships, or protein–protein interactions.

The app enables users to:

- Explore enriched pathway modules.
- Identify genes central to disease-associated biological processes.
- Compare pathway networks across increasing histological severity.
- Examine sex-aware pathway enrichment patterns.
- Evaluate protein-level connectivity using STRING.
- Browse individual gene-level meta-analysis results.
- Upload custom gene lists for comparison with MASLD meta-analysis results.

---

## Main Features

MASLD META NETWORK provides the following features:

- Interactive pathway–pathway network visualization.
- Pathway filtering by adjusted q-value.
- Edge pruning by Jaccard index or number of shared genes.
- Network clustering using Louvain, edge betweenness, or label propagation.
- User-defined naming of pathway clusters.
- Cluster–gene network visualization.
- Gene centrality analysis using degree, closeness, betweenness, or PageRank.
- Combined scoring based on meta-analysis score and network centrality.
- STRING protein–protein interaction network visualization.
- Stage-dependent pathway network comparison.
- Sex-aware pathway enrichment visualization.
- Gene-level transcriptome result browsing.
- Custom gene-list upload.
- Optional upload of logFC values.
- Ortholog conversion from zebrafish or mouse gene symbols to human gene symbols using DIOPT.
- Downloadable summary tables and interactive HTML networks.

---

## App Structure

The application is organized into six main tabs:

1. Top-Score Network
2. Process-Gene Network
3. STRING Network
4. Stage-Dependent Network
5. Sex-Aware Network
6. Transcriptome Browser

Each tab is designed to answer a different biological or analytical question.

---

## 1. Top-Score Network

The **Top-Score Network** tab visualizes pathway–pathway networks generated from over-representation analysis of top-scoring genes.

In this network, each node represents an enriched pathway, and each edge represents genes shared between two pathways. This structure allows users to identify groups of related biological processes and examine how enriched pathways are connected through shared gene content.

The sidebar allows users to select:

- The histological scoring system.
- Whether to include all genes or only genes with a specific direction of change.
- The MSigDB gene-set collection used for enrichment analysis.

Users can further filter enriched pathways according to q-value. Edges between pathways can be pruned using either the Jaccard index or the number of shared genes. The pruning threshold is adjustable by the user.

The tab supports pathway clustering using:

- Louvain clustering
- Edge betweenness clustering
- Label propagation clustering

After clustering, users can assign custom names to the resulting pathway clusters. Additional options allow users to control network movement and remove unconnected nodes.

In the network:

- Nodes represent pathways.
- Edges represent shared genes.
- Node size reflects enrichment significance.
- Edge color indicates the dominant direction of regulation.
- Green edges indicate predominantly upregulated genes.
- Red edges indicate predominantly downregulated genes.

Hovering over nodes and edges displays additional information, including pathway size, Jaccard index, number of shared genes, number of upregulated genes, and number of downregulated genes.

A summary table and an HTML version of the visualized network can be downloaded from the tab.

---

## 2. Process-Gene Network

The **Process-Gene Network** tab allows users to further investigate the pathway clusters identified in the Top-Score Network tab.

This tab creates a cluster–gene network in which pathway clusters are connected to the genes contributing to them. This allows users to identify genes that are central to specific biological modules rather than evaluating genes only as isolated meta-analysis hits.

The sidebar allows users to:

- Select which clusters to visualize.
- Apply additional q-value filtering.
- Select a gene centrality metric.
- Calculate a combined gene score.
- Adjust the relative weight of meta-analysis score and centrality score.
- Filter genes using a combined score threshold.

Available centrality metrics include:

- Degree
- Closeness
- Betweenness
- PageRank

The combined score integrates the gene’s total meta-analysis score with its centrality in the selected network. This helps prioritize genes that are both strongly supported by the meta-analysis and topologically important within a biological module.

Users can also upload a custom gene list as a TXT or CSV file. The uploaded file may contain only gene symbols or may include gene symbols with logFC values as a second column.

Gene symbols from zebrafish (*Danio rerio*) or mouse (*Mus musculus*) can be automatically converted to human orthologs using DIOPT.

When logFC values are provided, the app can filter genes according to concordance between the direction of change in the uploaded dataset and the direction observed in the MASLD meta-analysis.

Users can click on each cluster node to generate a subnetwork showing the genes and pathways involved in that specific cluster. The **Back to full network** button returns the user to the full cluster–gene network.

Hovering over gene nodes displays gene-level information, including:

- Gene name
- Total meta-analysis score
- Upregulated score component
- Downregulated score component
- Centrality value
- Combined score

The summary table and the interactive network can be downloaded.

---

## 3. STRING Network

The **STRING Network** tab visualizes protein–protein interaction networks derived from top-scoring genes identified in the MASLD meta-analysis.

This tab allows users to evaluate whether prioritized genes form connected protein-level modules and to identify genes that may be central within known interaction networks.

The sidebar allows users to select:

- The histological scoring system.
- The minimum STRING combined interaction score.
- The minimum total meta-analysis score required for displayed genes.
- The centrality metric used for network analysis.
- The relative weighting of centrality and meta-analysis score in the combined score.
- The combined score threshold used for filtering.

Available centrality metrics include:

- Degree
- Closeness
- Betweenness
- PageRank

Users can also upload custom gene lists with optional logFC values. Gene symbols from zebrafish or mouse can be converted to human orthologs using DIOPT.

If logFC values are provided, genes can be filtered according to concordance between the uploaded data and the MASLD meta-analysis direction.

The STRING Network tab also supports network clustering using:

- Louvain clustering
- Edge betweenness clustering
- Label propagation clustering

Genes can additionally be filtered according to their PubMed query status at the time of application release.

Advanced filters allow users to control the individual evidence channels contributing to the STRING combined score. This makes it possible to focus on specific types of protein interaction evidence.

A summary table and an HTML version of the visualized network are available for download.

---

## 4. Stage-Dependent Network

The **Stage-Dependent Network** tab visualizes how pathway networks change across increasing histological disease severity.

For each scoring system, the least severe group is used as the reference group. For NAS, the base group includes NAS 0, 1, and 2. For fibrosis, the base group includes F0 and F1.

Users can select a more advanced histological score group and visualize pathway networks derived from comparisons between that group and the base group.

This tab uses the same general network structure as the Top-Score Network tab:

- Nodes represent enriched pathways.
- Edges represent shared genes.
- Node size reflects enrichment significance.
- Edge color reflects the dominant direction of regulation.

The sidebar allows users to select:

- Histological scoring system.
- Direction of gene regulation.
- MSigDB gene-set collection.
- q-value threshold.
- Edge pruning method.
- Network clustering method.

A key feature of this tab is that pathways present in the previous disease stage but absent from the currently selected stage are shown as gray nodes. This helps users identify pathways that are gained, retained, or lost as disease severity increases.

This tab is designed to characterize how the enriched pathway space is remodeled during MASLD progression.

A summary table and an HTML version of the visualized network are available for download.

---

## 5. Sex-Aware Network

The **Sex-Aware Network** tab allows users to explore biological sex-associated differences in MASLD pathway enrichment.

For this analysis, the meta-analysis is repeated in patients for whom biological sex information is available. The resulting sex-stratified results are then filtered using the previously defined top-scoring genes.

The sidebar allows users to select:

- The histological scoring system.
- Direction of gene regulation.
- MSigDB gene-set collection.
- q-value threshold.
- Edge pruning method.
- Network clustering method.

In this network:

- Nodes represent pathways.
- Edges represent shared genes.
- Edge color indicates the dominant direction of regulation.
- Green edges indicate predominantly upregulated genes.
- Red edges indicate predominantly downregulated genes.
- Node color indicates whether the pathway is enriched in male patients, female patients, or both.
- Two-color nodes indicate pathways detected in both sex-stratified analyses.

Users can cluster the network using Louvain, edge betweenness, or label propagation algorithms. Resulting clusters can be named by the user.

Hovering over nodes and edges displays information such as pathway size, Jaccard index, number of shared genes, number of upregulated genes, and number of downregulated genes.

A summary table and an HTML version of the visualized network are available for download.

---

## 6. Transcriptome Browser

The **Transcriptome Browser** tab provides direct access to gene-level meta-analysis results.

Users can select a histological scoring system and a gene of interest from the sidebar. The app then displays general results for the selected gene across multiple analytical approaches.

The browser summarizes:

- Differential expression analysis of merged results.
- Individual dataset-level meta-analysis of pairwise score comparisons.
- Regression-based analyses.
- GeneDiseasePatterns analysis.
- Total gene score across methods.

Additional visualizations include:

- Dataset-specific regression results shown as boxplots.
- Pairwise differential expression results shown as a half-matrix of histological score comparisons.
- Progression-associated expression patterns generated using the GeneDiseasePatterns framework.

Users can also upload a TXT file containing a custom gene list. Gene symbols from zebrafish or mouse can be converted to human orthologs using DIOPT.

When a gene list is uploaded, the app generates a Gene List Dotplot. In this dotplot:

- Columns represent comparisons between the base group and increasing histological score groups.
- A dot indicates that the gene was significant in the corresponding individual differential expression analysis.
- Dot color represents the direction of regulation.
- Dot size represents the absolute median logFC value.

This tab allows users to validate individual genes identified in the network-based analyses.

---

## Recommended Tutorial Workflow

A typical analysis workflow is:

1. Open the **Top-Score Network** tab.
2. Select the histological scoring system of interest.
3. Select the gene direction and MSigDB gene-set collection.
4. Apply a q-value threshold.
5. Prune edges using Jaccard index or shared gene count.
6. Cluster the pathway network.
7. Assign biological names to the detected clusters.
8. Open the **Process-Gene Network** tab.
9. Select one or more pathway clusters.
10. Identify central genes using centrality and combined score.
11. Open the **STRING Network** tab.
12. Evaluate whether prioritized genes form connected protein-level modules.
13. Open the **Stage-Dependent Network** tab.
14. Examine whether pathway modules are gained, retained, or lost across disease progression.
15. Open the **Sex-Aware Network** tab.
16. Compare pathway enrichment patterns between male and female patients.
17. Use the **Transcriptome Browser** to inspect individual genes in detail.

---
