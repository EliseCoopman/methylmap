12/09/2022:
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
- def args_heatmap(): testing:
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