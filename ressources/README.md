# The resources folder

This folder contains mostly external data which are not novel but inherited from previous work.
The reference genome of the hamlet is from ___

Apart from this, the ressources folder also contains some intermediate results from this study that were reused a later stages. This refers specifically to the sub-folder images/ which contains the output folder from the image alignment software and 2 additional mask files that where created manually with a mean image of all aligned images in python and GIMP.

The following structure is assumed:

.
├── HP_genome_unmasked_01.dict
├── HP_genome_unmasked_01.fa
├── HP_genome_unmasked_01.fa.amb
├── HP_genome_unmasked_01.fa.ann
├── HP_genome_unmasked_01.fa.bwt
├── HP_genome_unmasked_01.fa.fai
├── HP_genome_unmasked_01.fa.gz
├── HP_genome_unmasked_01.fa.pac
├── HP_genome_unmasked_01.fa.sa
├── images
│   ├── body_mask.tif
│   ├── full_mask.tif
│   └── left_54off_59on
│       └── 3-registred
│           └── Modalities
│               └── RGB
│                   └── all
│                       ├── 28366gumboc-l1-s2-f1-c2-d1.png
│                       ├── 28377uniboc-l1-s1-f1-c2-d1.png
│                       ├── 28383uniboc-l1-s2-f1-c2-d0.png
│                       ├── 28384pueboc-l1-s2-f1-c2-d0.png
│                       ├── 28385nigboc-l1-s1-f1-c2-d0.png
│                       ├── 28386nigboc-l1-s1-f1-c2-d0.png
│                       ├── 28387nigboc-l1-s2-f1-c2-d0.png
│                       ├── 28388uniboc-l1-s1-f1-c2-d0.png
│                       ├── 28389abeboc-l1-s1-f1-c2-d0.png
│                       ├── 28390nigboc-l1-s2-f1-c2-d0.png
│                       ├── 28391uniboc-l1-s1-f1-c2-d0.png
│                       ├── 28392uniboc-l1-s2-f1-c2-d0.png
│                       ├── 28394nigboc-l1-s1-f1-c2-d0.png
│                       ├── 28399nigboc-l1-s1-f1-c2-d0.png
│                       ├── AG9RX46nigboc-l1-s2-f1-c2-d0.png
│                       ├── AG9RX48pueboc-l1-s1-f1-c2-d0.png
│                       ├── AG9RX49nigboc-l1-s1-f1-c2-d0.png
│                       ├── AG9RX50nigboc-l1-s2-f1-c2-d0.png
│                       ├── AG9RX51pueboc-l1-s1-f1-c2-d0.png
│                       ├── AG9RX53pueboc-l1-s2-f1-c2-d0.png
│                       ├── PL17_01uniboc-l1-s2-f1-c2-d0.png
│                       ├── PL17_02pueboc-l3-s1-f1-c2-d0.png
│                       ├── PL17_04pueboc-l1-s2-f1-c2-d0.png
│                       ├── PL17_05pueboc-l1-s2-f1-c2-d0.png
│                       ├── PL17_087nigbel-l4-s4-f4-c2-d1.png
│                       ├── PL17_088abebel-l3-s4-f4-c2-d1.png
│                       ├── PL17_089maybel-l3-s4-f4-c2-d1.png
│                       ├── PL17_090puebel-l1-s4-f4-c2-d1.png
│                       ├── PL17_091nigbel-l1-s4-f4-c2-d1.png
│                       ├── PL17_093nigbel-l1-s4-f4-c2-d1.png
│                       ├── PL17_094puebel-l1-s4-f4-c2-d1.png
│                       ├── PL17_095maybel-l11-s4-f4-c2-d1.png
│                       ├── PL17_096nigbel-l1-s4-f4-c2-d1.png
│                       ├── PL17_097indbel-l1-s4-f4-c2-d1.png
│                       ├── PL17_099indbel-l3-s4-f4-c2-d1.png
│                       ├── PL17_100indbel-l4-s4-f4-c2-d1.png
│                       ├── PL17_103puebel-l1-s4-f4-c2-d1.png
│                       ├── PL17_104nigbel-l1-s4-f4-c2-d1.png
│                       ├── PL17_105puebel-l3-s4-f4-c2-d1.png
│                       ├── PL17_106nigbel-l1-s4-f4-c2-d1.png
│                       ├── PL17_107puebel-l2-s4-f4-c2-d1.png
│                       ├── PL17_108nigbel-l1-s4-f4-c2-d1.png
│                       ├── PL17_109puebel-l2-s4-f4-c2-d1.png
│                       ├── PL17_110puebel-l1-s4-f4-c2-d1.png
│                       ├── PL17_111indbel-l1-s4-f4-c2-d1.png
│                       ├── PL17_112nigbel-l1-s4-f4-c2-d1.png
│                       ├── PL17_117puebel-l1-s4-f4-c2-d1.png
│                       ├── PL17_119maybel-l1-s4-f4-c2-d1.png
│                       ├── PL17_120maybel-l5-s4-f4-c2-d1.png
│                       ├── PL17_121maybel-l5-s4-f4-c2-d1.png
│                       ├── PL17_122maybel-l2-s4-f4-c2-d1.png
│                       ├── PL17_123maybel-l3-s4-f4-c2-d1.png
│                       ├── PL17_124maybel-l4-s4-f4-c2-d1.png
│                       ├── PL17_125tanbel-l7-s4-f4-c2-d1.png
│                       ├── PL17_126maybel-l1-s4-f4-c2-d1.png
│                       ├── PL17_127indbel-l6-s4-f4-c2-d1.png
│                       ├── PL17_128indbel-l1-s4-f4-c2-d1.png
│                       ├── PL17_132indbel-l2-s4-f4-c2-d1.png
│                       ├── PL17_134uniflo-l6-s4-f4-c2-d1.png
│                       ├── PL17_135uniflo-l6-s4-f4-c2-d1.png
│                       ├── PL17_136uniflo-l5-s4-f4-c2-d1.png
│                       ├── PL17_137uniflo-l4-s4-f4-c2-d1.png
│                       ├── PL17_138uniflo-r2-s4-f4-c2-d1.png
│                       ├── PL17_139pueflo-l1-s4-f4-c2-d1.png
│                       ├── PL17_140uniflo-l1-s4-f4-c2-d1.png
│                       ├── PL17_141uniflo-l8-s4-f4-c2-d1.png
│                       ├── PL17_142gemflo-l5-s4-f4-c2-d1.png
│                       ├── PL17_143uniflo-l4-s4-f4-c2-d1.png
│                       ├── PL17_144gemflo-l7-s4-f4-c2-d1.png
│                       ├── PL17_145gemflo-l10-s4-f4-c2-d1.png
│                       ├── PL17_148gemflo-l10-s4-f4-c2-d1.png
│                       ├── PL17_149nigflo-l1-s4-f4-c2-d1.png
│                       ├── PL17_153gemflo-l4-s4-f4-c2-d1.png
│                       ├── PL17_155pueflo-l1-s4-f4-c2-d1.png
│                       ├── PL17_157pueflo-l5-s4-f4-c2-d1.png
│                       ├── PL17_159pueflo-l5-s4-f4-c2-d1.png
│                       ├── PL17_160pueflo-l1-s4-f4-c2-d1.png
│                       ├── PL17_23nigpor-l1-s1-f3-c2-d1.png
│                       ├── PL17_35indpor-l1-s3-f3-c2-d1.png
│                       ├── PL17_37chlpor-l1-s3-f3-c2-d1.png
│                       ├── PL17_38chlpor-l1-s3-f3-c2-d1.png
│                       ├── PL17_39chlpor-l1-s3-f3-c2-d1.png
│                       ├── PL17_40chlpor-l1-s3-f3-c2-d1.png
│                       ├── PL17_41chlpor-l1-s3-f3-c2-d1.png
│                       ├── PL17_42chlpor-l1-s3-f3-c2-d1.png
│                       ├── PL17_43chlpor-l1-s3-f3-c2-d1.png
│                       ├── PL17_44chlpor-l1-s3-f3-c2-d1.png
│                       ├── PL17_50puepor-l2-s3-f3-c2-d1.png
│                       ├── PL17_53puepor-l1-s3-f3-c2-d1.png
│                       ├── PL17_54puepor-l2-s3-f3-c2-d1.png
│                       ├── PL17_55unipor-l1-s3-f3-c2-d1.png
│                       ├── PL17_56tanpor-l1-s3-f3-c2-d1.png
│                       ├── PL17_57puepor-l1-s3-f3-c2-d1.png
│                       ├── PL17_60puepor-l1-s3-f3-c2-d1.png
│                       ├── PL17_62puepor-l1-s3-f3-c2-d1.png
│                       ├── PL17_63unipor-l1-s3-f3-c2-d1.png
│                       ├── PL17_64indpor-l1-s3-f3-c2-d1.png
│                       ├── PL17_65puepor-l1-s3-f3-c2-d1.png
│                       ├── PL17_66unipor-l1-s3-f3-c2-d1.png
│                       ├── PL17_67unipor-l1-s3-f3-c2-d1.png
│                       ├── PL17_68gutpor-l1-s3-f3-c2-d1.png
│                       ├── PL17_69puepor-l1-s3-f3-c2-d1.png
│                       ├── PL17_70unipor-l1-s3-f3-c2-d1.png
│                       ├── PL17_71tanpor-l1-s3-f3-c2-d1.png
│                       ├── PL17_72tanpor-l1-s3-f3-c2-d1.png
│                       ├── PL17_73unipor-l1-s3-f3-c2-d1.png
│                       ├── PL17_74unipor-l1-s3-f3-c2-d1.png
│                       ├── PL17_75abepor-l1-s3-f3-c2-d1.png
│                       ├── PL17_76tanpor-l1-s3-f3-c2-d1.png
│                       ├── PL17_77unipor-l1-s3-f3-c2-d1.png
│                       ├── PL17_82puepor-l1-s3-f4-c2-d1.png
│                       ├── PL17_85indpor-l1-s3-f4-c2-d1.png
│                       └── PL17_86chlpor-l1-s3-f4-c2-d1.png
└── logos
    ├── belize.png
    ├── generic_hamlet_l.c.svg
    ├── generic_hamlet_r.c.svg
    ├── H_aberrans.l.cairo.png
    ├── H_aberrans.l.cairo.svg
    ├── H_atlahua.l.cairo.png
    ├── H_atlahua.l.cairo.svg
    ├── H_castroaguirrei.l.cairo.png
    ├── H_castroaguirrei.l.cairo.svg
    ├── H_chlorurus.l.cairo.png
    ├── H_chlorurus.l.cairo.svg
    ├── H_ecosur.l.cairo.png
    ├── H_ecosur.l.cairo.svg
    ├── H_floridae.l.cairo.png
    ├── H_floridae.l.cairo.svg
    ├── H_gemma.l.cairo.png
    ├── H_gemma.l.cairo.svg
    ├── H_gumigutta.l.cairo.png
    ├── H_gumigutta.l.cairo.svg
    ├── H_guttavarius.l.cairo.png
    ├── H_guttavarius.l.cairo.svg
    ├── H_indigo.l.cairo.png
    ├── H_indigo.l.cairo.svg
    ├── H_liberte.l.cairo.png
    ├── H_liberte.l.cairo.svg
    ├── H_maculiferus.l.cairo.svg
    ├── H_maculiferus.l.png
    ├── H_maya.l.cairo.png
    ├── H_maya.l.cairo.svg
    ├── H_nigricans.l.cairo.png
    ├── H_nigricans.l.cairo.svg
    ├── H_providencianus.l.cairo.png
    ├── H_providencianus.l.cairo.svg
    ├── H_puella.l.cairo.png
    ├── H_puella.l.cairo.svg
    ├── H_randallorum.l.cairo.png
    ├── H_randallorum.l.cairo.svg
    ├── H_sp.l.cairo.png
    ├── H_tan.l.cairo.png
    ├── H_tan.l.cairo.svg
    ├── H_unicolor.l.cairo.png
    ├── H_unicolor.l.cairo.svg
    ├── logo2.c.svg
    ├── logo.c.svg
    ├── pan.png
    ├── puer.png
    ├── S_tabacarius.l.cairo.svg
    ├── S_tabacarius.l.png
    ├── S_tigrinus.l.cairo.svg
    ├── S_tigrinus.l.png
    ├── S_tortugarum.l.cairo.svg
    ├── S_tortugarum.l.png
    └── us.png

7 directories, 177 files