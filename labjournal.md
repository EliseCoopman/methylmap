### V12/09/2022:
- --expand option
- def parse_annotation(): get information attributes: gene_name, transcript_id, locus_tag

13/09/2022:
- def parse_annotation(): problem when both locus_tag and gene_name present??
annotationfile["attributes"] = annotationfile.attributes.str.split(";")
    for x in annotationfile.attributes:
        for i in x:
            if i.startswith ("locus_tag") or i.startswith("gene_name"): #wat als beide aanwezig zijn?
                annotationfile['gene'] = i.split("=")[1]
            if i.startswith("transcript_id"):
                annotationfile['transcript'] = i.split("=")[1]
```python                
def args_heatmap(): testing:
    annot_axis = f'axis{annot_row}'
    num_methrows = len(table["MeanPosition"])
    annot_row = num_methrows
    if gtf:
        annotation_traces, y_max = gtf_annotation(gtf, region, simplify)
        for annot_trace in annotation_traces:
            heatmap.append_trace(trace=annot_trace, row=annot_row, col=1)
        heatmap["layout"][annot_axis].update(range=[-2,y_max + 1], showgrid=False,zeroline=False, showline=False, ticks='',showticklables=False)
        heatmap["layout"]["xaxis"].update(tickformat='g', seperatethousands=True,range=[region.begin,region.end])
        heatmap["layout"].update(barmode='overlay',title="Heatmap methylation frequency", hovermode='closest', plot_bgcolor='rgba(0,0,0,0)')
    html = heatmap.to_html()
```
20/09/2022:
- create subplots for annotation trace
- make heatmap using plotly graph_objects instead of plotly express to make subplots
- still a problem with displaying the calculation frequencies as text in the heatmap
- layout subplot (1,1) for the annotation track

21/09/2022:
- fix expand problem in annotation trace
- add def file_sniffer(), def is_gz_file(), def is_cram_file(), def is_bam_file(), def get_data(), def read_mods()
- write def parse_nanopolish(), rewrite args_heatmap() (not yet finished)
- problem: input overviewtable is one file, input nanopolish_calc_meth_freq or bam/cram are multiple files

