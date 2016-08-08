---
title: "Expression Report: Version 1"
author: "Shelby"
date: "08 August, 2016"
output:
  pdf_document:
    includes:
      in_header: header.tex
 
geometry: left=2cm, right = 2cm,  top = 2cm, bottom =2cm
---








#***Gene of Interest: GLS (ENSG00000115419)***

---

\begin{center}
\huge Expression Data
\end{center}

---


![](/home/shelby/for_andrew/ENSG00000115419SHH_Marked/ENSG00000115419_SHH_Marked_ExpressionReport_files/figure-latex/expr_subgroup-1.pdf)<!-- --> 


\newpage \begin{center}  \huge Expression vs. Clinical Features:  \newline \end{center}\small 1. Sex \newline\small 2. Age \newline\small 3. Stage \newline\small 4. Resection \newline\small 5. Status follow up \newline\small 6. Relapse \newline\small 7. M status at relapse \newline\small 8. Beta catenin IHC result \newline\small 9. Beta catenin Sequencing result \newline\small 10. TP53 IHC result \newline\small 11. TP53 Sequencing result \newline\small 12. MYC FISH result \newline\small 13. MYCN FISH result \newline\small 14. TERT Sequencing result \newline\small 15. TERT Taqman result \newline\small 16. Final pathology \newline\newpage \large Expression vs. Clinical Features:  \newline \newline \newline ![](/home/shelby/for_andrew/ENSG00000115419SHH_Marked/ENSG00000115419_SHH_Marked_ExpressionReport_files/figure-latex/clin_feat-1.pdf)<!-- --> 
![](/home/shelby/for_andrew/ENSG00000115419SHH_Marked/ENSG00000115419_SHH_Marked_ExpressionReport_files/figure-latex/clin_feat-2.pdf)<!-- --> 
![](/home/shelby/for_andrew/ENSG00000115419SHH_Marked/ENSG00000115419_SHH_Marked_ExpressionReport_files/figure-latex/clin_feat-3.pdf)<!-- --> 
![](/home/shelby/for_andrew/ENSG00000115419SHH_Marked/ENSG00000115419_SHH_Marked_ExpressionReport_files/figure-latex/clin_feat-4.pdf)<!-- --> 

\newpage

![](/home/shelby/for_andrew/ENSG00000115419SHH_Marked/ENSG00000115419_SHH_Marked_ExpressionReport_files/figure-latex/exprs_surv-1.pdf)<!-- --> 


\large Results from categorical and continuous COX modelling (all): \normalsize

                          n            HR  95% CI        P value
------------------------  -------  ------  -----------  --------
High Expression (Cat.)    68/222    1.028  0.638-1.66      0.911
High Expression (Cont.)   68/222    0.963  0.738-1.26      0.784



![](/home/shelby/for_andrew/ENSG00000115419SHH_Marked/ENSG00000115419_SHH_Marked_ExpressionReport_files/figure-latex/exprs_surv-2.pdf)<!-- --> 


\large Results from categorical and continuous COX modelling (WNT): \normalsize

                          n          HR  95% CI         P value
------------------------  -----  ------  ------------  --------
High Expression (Cat.)    2/27    0.742  0.0464-11.9      0.833
High Expression (Cont.)   2/27    0.674  0.013-34.8       0.845



![](/home/shelby/for_andrew/ENSG00000115419SHH_Marked/ENSG00000115419_SHH_Marked_ExpressionReport_files/figure-latex/exprs_surv-3.pdf)<!-- --> 


\large Results from categorical and continuous COX modelling (SHH): \normalsize

                          n           HR  95% CI        P value
------------------------  ------  ------  -----------  --------
High Expression (Cat.)    16/56    1.522  0.551-4.21      0.418
High Expression (Cont.)   16/56    1.558  0.769-3.16      0.218



![](/home/shelby/for_andrew/ENSG00000115419SHH_Marked/ENSG00000115419_SHH_Marked_ExpressionReport_files/figure-latex/exprs_surv-4.pdf)<!-- --> 


\large Results from categorical and continuous COX modelling (Grp3): \normalsize

                          n           HR  95% CI        P value
------------------------  ------  ------  -----------  --------
High Expression (Cat.)    25/53    1.118  0.509-2.46      0.781
High Expression (Cont.)   25/53    1.051  0.671-1.65      0.829



![](/home/shelby/for_andrew/ENSG00000115419SHH_Marked/ENSG00000115419_SHH_Marked_ExpressionReport_files/figure-latex/exprs_surv-5.pdf)<!-- --> 


\large Results from categorical and continuous COX modelling (Grp4): \normalsize

                          n           HR  95% CI        P value
------------------------  ------  ------  -----------  --------
High Expression (Cat.)    25/86    0.672  0.295-1.53      0.345
High Expression (Cont.)   25/86    0.688  0.379-1.25      0.219



\blandscape
![](/home/shelby/for_andrew/ENSG00000115419SHH_Marked/ENSG00000115419_SHH_Marked_ExpressionReport_files/figure-latex/isoform_plot-1.pdf)<!-- --> 
\elandscape

![](/home/shelby/for_andrew/ENSG00000115419SHH_Marked/ENSG00000115419_SHH_Marked_ExpressionReport_files/figure-latex/isoform_multi-1.pdf)<!-- --> ![](/home/shelby/for_andrew/ENSG00000115419SHH_Marked/ENSG00000115419_SHH_Marked_ExpressionReport_files/figure-latex/isoform_multi-2.pdf)<!-- --> ![](/home/shelby/for_andrew/ENSG00000115419SHH_Marked/ENSG00000115419_SHH_Marked_ExpressionReport_files/figure-latex/isoform_multi-3.pdf)<!-- --> ![](/home/shelby/for_andrew/ENSG00000115419SHH_Marked/ENSG00000115419_SHH_Marked_ExpressionReport_files/figure-latex/isoform_multi-4.pdf)<!-- --> 

\blandscape
![](/home/shelby/for_andrew/ENSG00000115419SHH_Marked/ENSG00000115419_SHH_Marked_ExpressionReport_files/figure-latex/isoform_bars-1.pdf)<!-- --> 
\elandscape

\blandscape

\elandscape

\blandscape

![](/home/shelby/for_andrew/ENSG00000115419SHH_Marked/ENSG00000115419_SHH_Marked_ExpressionReport_files/figure-latex/trans_by_iso-1.pdf)<!-- --> 
\elandscape



\begin{center}\huge Variants for: GLS\end{center}\normalsize

                                                                         Sample                                                                      Chromosome      Start         End        Variant   
----  --------------------------------------------------------------------------------------------------------------------------------------------  ------------  -----------  -----------  ------------
SYN    [NMB169](http://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&position=chr2%3A191743598-191747598&hgsid=447717037_yn6TgO7346K7iAggW1vSYaUVRegQ)       chr2       191745598    191745598    - : GCAGCA 
SYN    [NMB732](http://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&position=chr2%3A191796812-191800812&hgsid=447717037_yn6TgO7346K7iAggW1vSYaUVRegQ)       chr2       191798812    191798812      C : T    
SYN    [NMB803](http://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&position=chr2%3A191796812-191800812&hgsid=447717037_yn6TgO7346K7iAggW1vSYaUVRegQ)       chr2       191798812    191798812      C : T    
SYN    [NMB836](http://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&position=chr2%3A191796954-191800954&hgsid=447717037_yn6TgO7346K7iAggW1vSYaUVRegQ)       chr2       191798954    191798954      T : C    
SYN    [NMB433](http://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&position=chr2%3A191797592-191801592&hgsid=447717037_yn6TgO7346K7iAggW1vSYaUVRegQ)       chr2       191799592    191799592      A : G    
SYN    [NMB80](http://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&position=chr2%3A191797619-191801619&hgsid=447717037_yn6TgO7346K7iAggW1vSYaUVRegQ)        chr2       191799619    191799619      C : A    
SYN    [NMB439](http://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&position=chr2%3A191826326-191830326&hgsid=447717037_yn6TgO7346K7iAggW1vSYaUVRegQ)       chr2       191828326    191828326      G : C    



<!-- FROM HERE ARE THE TSS VS SUBGROUP --> 


\blandscape
![](/home/shelby/for_andrew/ENSG00000115419SHH_Marked/ENSG00000115419_SHH_Marked_ExpressionReport_files/figure-latex/tss_plot-1.pdf)<!-- --> 
\elandscape


![](/home/shelby/for_andrew/ENSG00000115419SHH_Marked/ENSG00000115419_SHH_Marked_ExpressionReport_files/figure-latex/tss_multi-1.pdf)<!-- --> ![](/home/shelby/for_andrew/ENSG00000115419SHH_Marked/ENSG00000115419_SHH_Marked_ExpressionReport_files/figure-latex/tss_multi-2.pdf)<!-- --> ![](/home/shelby/for_andrew/ENSG00000115419SHH_Marked/ENSG00000115419_SHH_Marked_ExpressionReport_files/figure-latex/tss_multi-3.pdf)<!-- --> 

\blandscape
![](/home/shelby/for_andrew/ENSG00000115419SHH_Marked/ENSG00000115419_SHH_Marked_ExpressionReport_files/figure-latex/tss_bars-1.pdf)<!-- --> 
\elandscape

\blandscape

![](/home/shelby/for_andrew/ENSG00000115419SHH_Marked/ENSG00000115419_SHH_Marked_ExpressionReport_files/figure-latex/trans_by_tss-1.pdf)<!-- --> 
\elandscape



<!-- FROM HERE ARE THE CDS VS SUBGROUP --> 


\blandscape
![](/home/shelby/for_andrew/ENSG00000115419SHH_Marked/ENSG00000115419_SHH_Marked_ExpressionReport_files/figure-latex/cds_plot-1.pdf)<!-- --> 
\elandscape

![](/home/shelby/for_andrew/ENSG00000115419SHH_Marked/ENSG00000115419_SHH_Marked_ExpressionReport_files/figure-latex/cds_multi-1.pdf)<!-- --> ![](/home/shelby/for_andrew/ENSG00000115419SHH_Marked/ENSG00000115419_SHH_Marked_ExpressionReport_files/figure-latex/cds_multi-2.pdf)<!-- --> 
![](/home/shelby/for_andrew/ENSG00000115419SHH_Marked/ENSG00000115419_SHH_Marked_ExpressionReport_files/figure-latex/cds_multi-3.pdf)<!-- --> 

\blandscape
![](/home/shelby/for_andrew/ENSG00000115419SHH_Marked/ENSG00000115419_SHH_Marked_ExpressionReport_files/figure-latex/cds_bars-1.pdf)<!-- --> 
\elandscape

\blandscape
![](/home/shelby/for_andrew/ENSG00000115419SHH_Marked/ENSG00000115419_SHH_Marked_ExpressionReport_files/figure-latex/trans_by_cds-1.pdf)<!-- --> 
\elandscape

\blandscape


\begin{center}\huge Fusions for: GLS\end{center}\footnotesize\large No Fusions Found


\elandscape




