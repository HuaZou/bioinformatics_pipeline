rule summarize_metadata:
    input:
        metadata = config["metadata"]
    output:
        os.path.join(config["result"]["report"], "metadata_summary.md")
    run:
        df = pd.read_csv(input.metadata, sep='\t')
        cols = df.columns
        df2 = pd.DataFrame(columns =[0,1], index=cols)
        for col in cols:
            if col in df.columns:
                vc = df[col].value_counts()
                if vc.index.shape == (0,):
                    df2.loc[col, 0] = '(no values in column)'
                    df2.loc[col, 1] = '--'
                else:
                    df2.loc[col, 0] = vc.index[0]
                    df2.loc[col, 1] = vc.values[0]
            else:
                df2.loc[col, 0] = '(column not provided)'
                df2.loc[col, 1] = '--'
        df2.columns = ['Most common value', 'Count']
        df2.index.name = 'Column name'
        outstr = tabulate(df2, tablefmt="pipe", headers="keys")
        with open(output[0], 'w') as target:
            target.write(outstr)
            target.write('\n')

rule generate_report_md:
    input:
        configfile   = "config.yaml",
        mdsummary    = os.path.join(config["result"]["report"], "metadata_summary.md"),
        fastq        = os.path.join(config["result"]["remove"], "paired-end-demux-trim-fastq_counts_describe.md"), 
        amplicontype = os.path.join(config["result"]["denoise"], "repseqs_amplicon_type.txt"),
        repseqstsv   = os.path.join(config["result"]["tree"], "repseqs_properties.tsv"),
        repseqspdf   = os.path.join(config["result"]["tree"], "repseqs_properties.pdf"),
        repseqstree  = os.path.join(config["result"]["tree"], "rooted_tree.qzv"),
        repseqsoutliers   = os.path.join(config["result"]["tree"], "repseqs_to_filter_outliers.tsv"),
        repseqsunassigned = os.path.join(config["result"]["tree"], "repseqs_to_filter_unassigned.tsv"),
        samples      = os.path.join(config["result"]["denoise"], "table_summary_samples.txt"),
        features     = os.path.join(config["result"]["denoise"], "table_summary_features.txt"),
        visfastq     = os.path.join(config["result"]["remove"], "paired-end-demux-trim-fastq_summary.qzv"),
        visrepseqs   = os.path.join(config["result"]["denoise"], "repseqs_final.qzv"),
        tsvtaxonomy  = os.path.join(config["result"]["taxonomy"], "taxonomy.tsv"), 
        vistaxonomy  = os.path.join(config["result"]["taxonomy"], "taxonomy.qzv"),
        vistable     = os.path.join(config["result"]["denoise"], "table_final.qzv"),
        vistaxbar    = os.path.join(config["result"]["taxonomy"], "taxa_barplot.qzv"),
        visalpharare = os.path.join(config["result"]["alpha"], "alpha_rarefaction.qzv"),
        visevengs    = os.path.join(config["result"]["alpha"], "evenness_group_significance.qzv"),
        visfaithgs   = os.path.join(config["result"]["alpha"], "faith_pd_group_significance.qzv"),
        visobsfeaturesgs = os.path.join(config["result"]["alpha"], "observed_features_group_significance.qzv"),
        visshannongs = os.path.join(config["result"]["alpha"], "shannon_group_significance.qzv"),
        visbcemp     = os.path.join(config["result"]["beta"], "bray_curtis_emperor.qzv"),
        visbcgs      = os.path.join(config["result"]["beta"], "bray_curtis_group_significance.qzv"),
        visjacemp    = os.path.join(config["result"]["beta"], "jaccard_emperor.qzv"),
        visjacgs     = os.path.join(config["result"]["beta"], "jaccard_group_significance.qzv"),
        viswuemp     = os.path.join(config["result"]["beta"], "weighted_unifrac_emperor.qzv"),
        viswugs      = os.path.join(config["result"]["beta"], "weighted_unifrac_group_significance.qzv"),
        visuwuemp    = os.path.join(config["result"]["beta"], "unweighted_unifrac_emperor.qzv"),
        visuwugs     = os.path.join(config["result"]["beta"], "unweighted_unifrac_group_significance.qzv")
    log:
        os.path.join(config["logs"]["report"], "generate_report_markdown.log")     
    output:
        os.path.join(config["result"]["report"], "report.md")
    shell:
        "echo '# Report' > {output};"
        "echo '' >> {output};"
        "echo 'View this HTML report with [Chrome](https://www.google.com/chrome/){{target=\"_blank\"}} for best results.' >> {output};"
        "echo '' >> {output};"
        "echo 'To view the linked files below: ' >> {output};"
        "echo '' >> {output};"
        "echo '* QZV (QIIME 2 visualization): click to download, then drag and drop in [https://view.qiime2.org](https://view.qiime2.org){{target=\"_blank\"}}.' >> {output};"
        "echo '* TSV (tab-separated values): click to download, then open in Microsoft Excel or Tabview (command line tool that comes with Tourmaline).' >> {output};"
        "echo '* PDF (portable document format): click to open and view in new tab.' >> {output};"
        "echo '* Markdown and text: click to open and view in new tab.' >> {output};"
        "echo '' >> {output};"
        "echo 'Note: Downloaded files can be deleted after viewing because they are already stored in your Tourmaline directory.' >> {output};"
        "echo '' >> {output};"
        "echo '---' >> {output};"
        "echo '' >> {output};"
        "echo '## Metadata Summary' >> {output};"
        "echo '' >> {output};"
        "echo Markdown: \[{input.mdsummary}\]\(../{input.mdsummary}\){{target=\"_blank\"}} >> {output};"
        "echo '' >> {output};"
        "cat {input.mdsummary} >> {output};"
        "echo '' >> {output};"
        "echo '' >> {output};"
        "echo '---' >> {output};"
        "echo '' >> {output};"
        "echo '## Fastq Sequences Information' >> {output};"
        "echo '' >> {output};"
        "echo '### Fastq Sequences per Sample' >> {output};"
        "echo '' >> {output};"
        "echo Markdown: \[{input.fastq}\]\(../{input.fastq}\){{target=\"_blank\"}} >> {output};"
        "echo '' >> {output};"
        "cat {input.fastq} >> {output};"
        "echo '' >> {output};"
        "echo '### Visualization of Fastq Sequences' >> {output};"
        "echo '' >> {output};"
        "echo QZV: \[{input.visfastq}\]\(../{input.visfastq}\){{target=\"_blank\"}} >> {output};"
        "echo '' >> {output};"
        "echo '---' >> {output};"
        "echo '' >> {output};"
        "echo '## Representative Sequences Information' >> {output};"
        "echo '' >> {output};"
        "echo '### Representative Sequences Properties Table' >> {output};"
        "echo '' >> {output};"
        "echo TSV: \[{input.repseqstsv}\]\(../{input.repseqstsv}\){{target=\"_blank\"}} >> {output};"
        "echo '' >> {output};"
        "echo 'Columns:' >> {output};"
        "echo '' >> {output};"
        "echo '* featureid' >> {output};"
        "echo '* length - length (bp) not including gaps' >> {output};"
        "echo '* gaps - gaps (bp) in masked multiple sequence alignment' >> {output};"
        "echo '* outlier - outlier (True/False) determined by OD-seq' >> {output};"
        "echo '* taxonomy - domain level' >> {output};"
        "echo '* observations - total observations summed across all samples (unrarefied)' >> {output};"
        "echo '* log10(observations) - log base-10 of total observations' >> {output};"
        "echo '' >> {output};"
        "echo '### Representative Sequences Properties Plot' >> {output};"
        "echo '' >> {output};"
        "echo PDF: \[{input.repseqspdf}\]\(../{input.repseqspdf}\){{target=\"_blank\"}} >> {output};"
        "echo '' >> {output};"
        "echo 'Plot elements:' >> {output};"
        "echo '' >> {output};"
        "echo '* x: length (bp) not including gaps' >> {output};"
        "echo '* y: gaps (bp) in masked multiple sequence alignment' >> {output};"
        "echo '* color: taxonomy (domain)' >> {output};"
        "echo '* size: log10(observations)' >> {output};"
        "echo '* facets: outlier (True/False) determined by OD-seq' >> {output};"
        "echo '' >> {output};"
        "echo '### Representative Sequences Rooted Tree' >> {output};"
        "echo '' >> {output};"
        "echo QZV: \[{input.repseqstree}\]\(../{input.repseqstree}\){{target=\"_blank\"}} >> {output};"
        "echo '' >> {output};"
        "echo '### Representative Sequences to Filter' >> {output};"
        "echo '' >> {output};"
        "echo 'To filter, go to 02-output-{{method}}-{{filter}}/02-alignment-tree, merge or copy repseqs_to_filter_outliers.tsv and/or repseqs_to_filter_unassigned.tsv, rename to 00-data/repseqs_to_filter_{{method}}.tsv, then run Snakemake in filtered mode.' >> {output};"
        "echo '' >> {output};"
        "echo Outliers TSV: \[{input.repseqsoutliers}\]\(../{input.repseqsoutliers}\){{target=\"_blank\"}} >> {output};"
        "echo '' >> {output};"
        "echo Unassigned TSV: \[{input.repseqsunassigned}\]\(../{input.repseqsunassigned}\){{target=\"_blank\"}} >> {output};"
        "echo '' >> {output};"
        "echo '### Visualization of Representative Sequences' >> {output};"
        "echo '' >> {output};"
        "echo QZV: \[{input.visrepseqs}\]\(../{input.visrepseqs}\){{target=\"_blank\"}} >> {output};"
        "echo '' >> {output};"
        "echo '### Taxonomy Table' >> {output};"
        "echo '' >> {output};"
        "echo TSV: \[{input.tsvtaxonomy}\]\(../{input.tsvtaxonomy}\){{target=\"_blank\"}} >> {output};"
        "echo '' >> {output};"
        "echo '### Visualization of Taxonomy' >> {output};"
        "echo '' >> {output};"
        "echo QZV: \[{input.vistaxonomy}\]\(../{input.vistaxonomy}\){{target=\"_blank\"}} >> {output};"
        "echo '' >> {output};"
        "echo '### Predicted Amplicon Type' >> {output};"
        "echo '' >> {output};"
        "echo Text: \[{input.amplicontype}\]\(../{input.amplicontype}\){{target=\"_blank\"}} >> {output};"
        "echo '' >> {output};"
        "echo '```' >> {output};"
        "head -n 5 {input.amplicontype} >> {output};"
        "echo '```' >> {output};"
        "echo '' >> {output};"
        "echo '---' >> {output};"
        "echo '' >> {output};"
        "echo '## Observation Table Information' >> {output};"
        "echo '' >> {output};"
        "echo '### Table Summary: Samples' >> {output};"
        "echo '' >> {output};"
        "echo Text: \[{input.samples}\]\(../{input.samples}\){{target=\"_blank\"}} >> {output};"
        "echo '' >> {output};"
        "echo '```' >> {output};"
        "head -n 13 {input.samples} >> {output};"
        "echo '```' >> {output};"
        "echo '' >> {output};"
        "echo '### Table Summary: Features' >> {output};"
        "echo '' >> {output};"
        "echo Text: \[{input.features}\]\(../{input.features}\){{target=\"_blank\"}} >> {output};"
        "echo '' >> {output};"
        "echo '```' >> {output};"
        "head -n 13 {input.features} >> {output};"
        "echo '```' >> {output};"
        "echo '' >> {output};"
        "echo '### Visualization of Table Summary' >> {output};"
        "echo '' >> {output};"
        "echo QZV: \[{input.vistable}\]\(../{input.vistable}\){{target=\"_blank\"}} >> {output};"
        "echo '' >> {output};"
        "echo '---' >> {output};"
        "echo '' >> {output};"
        "echo '## Taxonomic Diversity Results' >> {output};"
        "echo '' >> {output};"
        "echo '### Taxonomy Barplot' >> {output};"
        "echo '' >> {output};"
        "echo QZV: \[{input.vistaxbar}\]\(../{input.vistaxbar}\){{target=\"_blank\"}} >> {output};"
        "echo '' >> {output};"
        "echo '---' >> {output};"
        "echo '' >> {output};"
        "echo '## Alpha Diversity Results' >> {output};"
        "echo '' >> {output};"
        "echo '### Alpha Rarefaction' >> {output};"
        "echo '' >> {output};"
        "echo QZV: \[{input.visalpharare}\]\(../{input.visalpharare}\){{target=\"_blank\"}} >> {output};"
        "echo '' >> {output};"
        "echo '### Alpha Group Significance' >> {output};"
        "echo '' >> {output};"
        "echo Evenness QZV: \[{input.visevengs}\]\(../{input.visevengs}\){{target=\"_blank\"}} >> {output};"
        "echo '' >> {output};"
        "echo Faith PD QZV: \[{input.visfaithgs}\]\(../{input.visfaithgs}\){{target=\"_blank\"}} >> {output};"
        "echo '' >> {output};"
        "echo Observed features QZV: \[{input.visobsfeaturesgs}\]\(../{input.visobsfeaturesgs}\){{target=\"_blank\"}} >> {output};"
        "echo '' >> {output};"
        "echo Shannon QZV: \[{input.visshannongs}\]\(../{input.visshannongs}\){{target=\"_blank\"}} >> {output};"
        "echo '' >> {output};"
        "echo '---' >> {output};"
        "echo '' >> {output};"
        "echo '## Beta Diversity Results' >> {output};"
        "echo '' >> {output};"
        "echo '### PCoA Emperor Plots' >> {output};"
        "echo '' >> {output};"
        "echo Bray-Curtis QZV: \[{input.visbcemp}\]\(../{input.visbcemp}\){{target=\"_blank\"}} >> {output};"
        "echo '' >> {output};"
        "echo Jaccard QZV: \[{input.visjacemp}\]\(../{input.visjacemp}\){{target=\"_blank\"}} >> {output};"
        "echo '' >> {output};"
        "echo Weighted UniFrac QZV: \[{input.viswuemp}\]\(../{input.viswuemp}\){{target=\"_blank\"}} >> {output};"
        "echo '' >> {output};"
        "echo Unweighted UniFrac QZV: \[{input.visuwuemp}\]\(../{input.visuwuemp}\){{target=\"_blank\"}} >> {output};"
        "echo '' >> {output};"
        "echo '### Beta Group Significance' >> {output};"
        "echo '' >> {output};"
        "echo Bray-Curtis QZV: \[{input.visbcgs}\]\(../{input.visbcgs}\){{target=\"_blank\"}} >> {output};"
        "echo '' >> {output};"
        "echo Jaccard QZV: \[{input.visjacgs}\]\(../{input.visjacgs}\){{target=\"_blank\"}} >> {output};"
        "echo '' >> {output};"
        "echo Weighted UniFrac QZV: \[{input.viswugs}\]\(../{input.viswugs}\){{target=\"_blank\"}} >> {output};"
        "echo '' >> {output};"
        "echo Unweighted UniFrac QZV: \[{input.visuwugs}\]\(../{input.visuwugs}\){{target=\"_blank\"}} >> {output};"
        "echo '' >> {output};"
        "echo '---' >> {output};"
        "echo '' >> {output};"
        "echo '## Tourmaline Config File' >> {output};"
        "echo '' >> {output};"
        "echo YAML: \[{input.configfile}\]\(../{input.configfile}\){{target=\"_blank\"}} >> {output};"
        "echo '' >> {output};"
        "echo '```' >> {output};"
        "cat {input.configfile} >> {output};"
        "echo '```' >> {output}; 2>{log}"

rule generate_report_html:
    input:
        os.path.join(config["result"]["report"], "report.md")
    params:
        theme = config["params"]["report"]["report_theme"]
    log:
        os.path.join(config["logs"]["report"], "generate_report_html.log") 
    output:
        os.path.join(config["result"]["report"], "report.html")
    shell:
        "pandoc -i {input} -o {output};"
        "echo '<!DOCTYPE html>' > header.html;"
        "echo '<html>' >> header.html;"
        "echo '<head>' >> header.html;"
        "echo '<link rel=\"stylesheet\" type=\"text/css\" href=\"../css/{params.theme}.css\">' >> header.html;"
        "echo '</head>' >> header.html;"
        "echo '<body>' >> header.html;"
        "echo '' >> header.html;"
        "cat header.html {output} > temp.html;"
        "echo '' >> temp.html;"
        "echo '</body>' >> temp.html;"
        "echo '</html>' >> temp.html;"
        "mv temp.html {output};"
        "rm header.html 2>{log}"
