# The resources folder

This folder contains mostly external data which are not novel but inherited from previous work.
The reference genome of the hamlet is from ___

Apart from this, the ressources folder also contains some intermediate results from this study that were reused a later stages. This refers specifically to the sub-folder images/ which contains the output folder from the image alignment software and 2 additional mask files that where created manually with a mean image of all aligned images in python and GIMP.

The following structure is assumed:

```
.<br>
├── HP_genome_unmasked_01.dict<br>
├── HP_genome_unmasked_01.fa<br>
├── HP_genome_unmasked_01.fa.amb<br>
├── HP_genome_unmasked_01.fa.ann<br>
├── HP_genome_unmasked_01.fa.bwt<br>
├── HP_genome_unmasked_01.fa.fai<br>
├── HP_genome_unmasked_01.fa.gz<br>
├── HP_genome_unmasked_01.fa.pac<br>
├── HP_genome_unmasked_01.fa.sa<br>
├── images<br>
│   ├── body_mask.tif<br>
│   ├── full_mask.tif<br>
│   └── left_54off_59on<br>
│       └── 3-registred<br>
│           └── Modalities<br>
│               └── RGB<br>
│                   └── all<br>
│                       ├── 28366gumboc-l1-s2-f1-c2-d1.png<br>
│                       ├── 28377uniboc-l1-s1-f1-c2-d1.png<br>
│                       ├── 28383uniboc-l1-s2-f1-c2-d0.png<br>
│                       ├── 28384pueboc-l1-s2-f1-c2-d0.png<br>
│                       ├── 28385nigboc-l1-s1-f1-c2-d0.png<br>
│                       ├── 28386nigboc-l1-s1-f1-c2-d0.png<br>
│                       ├── 28387nigboc-l1-s2-f1-c2-d0.png<br>
│                       ├── 28388uniboc-l1-s1-f1-c2-d0.png<br>
│                       ├── 28389abeboc-l1-s1-f1-c2-d0.png<br>
│                       ├── 28390nigboc-l1-s2-f1-c2-d0.png<br>
│                       ├── 28391uniboc-l1-s1-f1-c2-d0.png<br>
│                       ├── 28392uniboc-l1-s2-f1-c2-d0.png<br>
│                       ├── 28394nigboc-l1-s1-f1-c2-d0.png<br>
│                       ├── 28399nigboc-l1-s1-f1-c2-d0.png<br>
│                       ├── AG9RX46nigboc-l1-s2-f1-c2-d0.png<br>
│                       ├── AG9RX48pueboc-l1-s1-f1-c2-d0.png<br>
│                       ├── AG9RX49nigboc-l1-s1-f1-c2-d0.png<br>
│                       ├── AG9RX50nigboc-l1-s2-f1-c2-d0.png<br>
│                       ├── AG9RX51pueboc-l1-s1-f1-c2-d0.png<br>
│                       ├── AG9RX53pueboc-l1-s2-f1-c2-d0.png<br>
│                       ├── PL17_01uniboc-l1-s2-f1-c2-d0.png<br>
│                       ├── PL17_02pueboc-l3-s1-f1-c2-d0.png<br>
│                       ├── PL17_04pueboc-l1-s2-f1-c2-d0.png<br>
│                       ├── PL17_05pueboc-l1-s2-f1-c2-d0.png<br>
│                       ├── PL17_087nigbel-l4-s4-f4-c2-d1.png<br>
│                       ├── PL17_088abebel-l3-s4-f4-c2-d1.png<br>
│                       ├── PL17_089maybel-l3-s4-f4-c2-d1.png<br>
│                       ├── PL17_090puebel-l1-s4-f4-c2-d1.png<br>
│                       ├── PL17_091nigbel-l1-s4-f4-c2-d1.png<br>
│                       ├── PL17_093nigbel-l1-s4-f4-c2-d1.png<br>
│                       ├── PL17_094puebel-l1-s4-f4-c2-d1.png<br>
│                       ├── PL17_095maybel-l11-s4-f4-c2-d1.png<br>
│                       ├── PL17_096nigbel-l1-s4-f4-c2-d1.png<br>
│                       ├── PL17_097indbel-l1-s4-f4-c2-d1.png<br>
│                       ├── PL17_099indbel-l3-s4-f4-c2-d1.png<br>
│                       ├── PL17_100indbel-l4-s4-f4-c2-d1.png<br>
│                       ├── PL17_103puebel-l1-s4-f4-c2-d1.png<br>
│                       ├── PL17_104nigbel-l1-s4-f4-c2-d1.png<br>
│                       ├── PL17_105puebel-l3-s4-f4-c2-d1.png<br>
│                       ├── PL17_106nigbel-l1-s4-f4-c2-d1.png<br>
│                       ├── PL17_107puebel-l2-s4-f4-c2-d1.png<br>
│                       ├── PL17_108nigbel-l1-s4-f4-c2-d1.png<br>
│                       ├── PL17_109puebel-l2-s4-f4-c2-d1.png<br>
│                       ├── PL17_110puebel-l1-s4-f4-c2-d1.png<br>
│                       ├── PL17_111indbel-l1-s4-f4-c2-d1.png<br>
│                       ├── PL17_112nigbel-l1-s4-f4-c2-d1.png<br>
│                       ├── PL17_117puebel-l1-s4-f4-c2-d1.png<br>
│                       ├── PL17_119maybel-l1-s4-f4-c2-d1.png<br>
│                       ├── PL17_120maybel-l5-s4-f4-c2-d1.png<br>
│                       ├── PL17_121maybel-l5-s4-f4-c2-d1.png<br>
│                       ├── PL17_122maybel-l2-s4-f4-c2-d1.png<br>
│                       ├── PL17_123maybel-l3-s4-f4-c2-d1.png<br>
│                       ├── PL17_124maybel-l4-s4-f4-c2-d1.png<br>
│                       ├── PL17_125tanbel-l7-s4-f4-c2-d1.png<br>
│                       ├── PL17_126maybel-l1-s4-f4-c2-d1.png<br>
│                       ├── PL17_127indbel-l6-s4-f4-c2-d1.png<br>
│                       ├── PL17_128indbel-l1-s4-f4-c2-d1.png<br>
│                       ├── PL17_132indbel-l2-s4-f4-c2-d1.png<br>
│                       ├── PL17_134uniflo-l6-s4-f4-c2-d1.png<br>
│                       ├── PL17_135uniflo-l6-s4-f4-c2-d1.png<br>
│                       ├── PL17_136uniflo-l5-s4-f4-c2-d1.png<br>
│                       ├── PL17_137uniflo-l4-s4-f4-c2-d1.png<br>
│                       ├── PL17_138uniflo-r2-s4-f4-c2-d1.png<br>
│                       ├── PL17_139pueflo-l1-s4-f4-c2-d1.png<br>
│                       ├── PL17_140uniflo-l1-s4-f4-c2-d1.png<br>
│                       ├── PL17_141uniflo-l8-s4-f4-c2-d1.png<br>
│                       ├── PL17_142gemflo-l5-s4-f4-c2-d1.png<br>
│                       ├── PL17_143uniflo-l4-s4-f4-c2-d1.png<br>
│                       ├── PL17_144gemflo-l7-s4-f4-c2-d1.png<br>
│                       ├── PL17_145gemflo-l10-s4-f4-c2-d1.png<br>
│                       ├── PL17_148gemflo-l10-s4-f4-c2-d1.png<br>
│                       ├── PL17_149nigflo-l1-s4-f4-c2-d1.png<br>
│                       ├── PL17_153gemflo-l4-s4-f4-c2-d1.png<br>
│                       ├── PL17_155pueflo-l1-s4-f4-c2-d1.png<br>
│                       ├── PL17_157pueflo-l5-s4-f4-c2-d1.png<br>
│                       ├── PL17_159pueflo-l5-s4-f4-c2-d1.png<br>
│                       ├── PL17_160pueflo-l1-s4-f4-c2-d1.png<br>
│                       ├── PL17_23nigpor-l1-s1-f3-c2-d1.png<br>
│                       ├── PL17_35indpor-l1-s3-f3-c2-d1.png<br>
│                       ├── PL17_37chlpor-l1-s3-f3-c2-d1.png<br>
│                       ├── PL17_38chlpor-l1-s3-f3-c2-d1.png<br>
│                       ├── PL17_39chlpor-l1-s3-f3-c2-d1.png<br>
│                       ├── PL17_40chlpor-l1-s3-f3-c2-d1.png<br>
│                       ├── PL17_41chlpor-l1-s3-f3-c2-d1.png<br>
│                       ├── PL17_42chlpor-l1-s3-f3-c2-d1.png<br>
│                       ├── PL17_43chlpor-l1-s3-f3-c2-d1.png<br>
│                       ├── PL17_44chlpor-l1-s3-f3-c2-d1.png<br>
│                       ├── PL17_50puepor-l2-s3-f3-c2-d1.png<br>
│                       ├── PL17_53puepor-l1-s3-f3-c2-d1.png<br>
│                       ├── PL17_54puepor-l2-s3-f3-c2-d1.png<br>
│                       ├── PL17_55unipor-l1-s3-f3-c2-d1.png<br>
│                       ├── PL17_56tanpor-l1-s3-f3-c2-d1.png<br>
│                       ├── PL17_57puepor-l1-s3-f3-c2-d1.png<br>
│                       ├── PL17_60puepor-l1-s3-f3-c2-d1.png<br>
│                       ├── PL17_62puepor-l1-s3-f3-c2-d1.png<br>
│                       ├── PL17_63unipor-l1-s3-f3-c2-d1.png<br>
│                       ├── PL17_64indpor-l1-s3-f3-c2-d1.png<br>
│                       ├── PL17_65puepor-l1-s3-f3-c2-d1.png<br>
│                       ├── PL17_66unipor-l1-s3-f3-c2-d1.png<br>
│                       ├── PL17_67unipor-l1-s3-f3-c2-d1.png<br>
│                       ├── PL17_68gutpor-l1-s3-f3-c2-d1.png<br>
│                       ├── PL17_69puepor-l1-s3-f3-c2-d1.png<br>
│                       ├── PL17_70unipor-l1-s3-f3-c2-d1.png<br>
│                       ├── PL17_71tanpor-l1-s3-f3-c2-d1.png<br>
│                       ├── PL17_72tanpor-l1-s3-f3-c2-d1.png<br>
│                       ├── PL17_73unipor-l1-s3-f3-c2-d1.png<br>
│                       ├── PL17_74unipor-l1-s3-f3-c2-d1.png<br>
│                       ├── PL17_75abepor-l1-s3-f3-c2-d1.png<br>
│                       ├── PL17_76tanpor-l1-s3-f3-c2-d1.png<br>
│                       ├── PL17_77unipor-l1-s3-f3-c2-d1.png<br>
│                       ├── PL17_82puepor-l1-s3-f4-c2-d1.png<br>
│                       ├── PL17_85indpor-l1-s3-f4-c2-d1.png<br>
│                       └── PL17_86chlpor-l1-s3-f4-c2-d1.png<br>
└── logos<br>
    ├── belize.png<br>
    ├── generic_hamlet_l.c.svg<br>
    ├── generic_hamlet_r.c.svg<br>
    ├── H_aberrans.l.cairo.png<br>
    ├── H_aberrans.l.cairo.svg<br>
    ├── H_atlahua.l.cairo.png<br>
    ├── H_atlahua.l.cairo.svg<br>
    ├── H_castroaguirrei.l.cairo.png<br>
    ├── H_castroaguirrei.l.cairo.svg<br>
    ├── H_chlorurus.l.cairo.png<br>
    ├── H_chlorurus.l.cairo.svg<br>
    ├── H_ecosur.l.cairo.png<br>
    ├── H_ecosur.l.cairo.svg<br>
    ├── H_floridae.l.cairo.png<br>
    ├── H_floridae.l.cairo.svg<br>
    ├── H_gemma.l.cairo.png<br>
    ├── H_gemma.l.cairo.svg<br>
    ├── H_gumigutta.l.cairo.png<br>
    ├── H_gumigutta.l.cairo.svg<br>
    ├── H_guttavarius.l.cairo.png<br>
    ├── H_guttavarius.l.cairo.svg<br>
    ├── H_indigo.l.cairo.png<br>
    ├── H_indigo.l.cairo.svg<br>
    ├── H_liberte.l.cairo.png<br>
    ├── H_liberte.l.cairo.svg<br>
    ├── H_maculiferus.l.cairo.svg<br>
    ├── H_maculiferus.l.png<br>
    ├── H_maya.l.cairo.png<br>
    ├── H_maya.l.cairo.svg<br>
    ├── H_nigricans.l.cairo.png<br>
    ├── H_nigricans.l.cairo.svg<br>
    ├── H_providencianus.l.cairo.png<br>
    ├── H_providencianus.l.cairo.svg<br>
    ├── H_puella.l.cairo.png<br>
    ├── H_puella.l.cairo.svg<br>
    ├── H_randallorum.l.cairo.png<br>
    ├── H_randallorum.l.cairo.svg<br>
    ├── H_sp.l.cairo.png<br>
    ├── H_tan.l.cairo.png<br>
    ├── H_tan.l.cairo.svg<br>
    ├── H_unicolor.l.cairo.png<br>
    ├── H_unicolor.l.cairo.svg<br>
    ├── logo2.c.svg<br>
    ├── logo.c.svg<br>
    ├── pan.png<br>
    ├── puer.png<br>
    ├── S_tabacarius.l.cairo.svg<br>
    ├── S_tabacarius.l.png<br>
    ├── S_tigrinus.l.cairo.svg<br>
    ├── S_tigrinus.l.png<br>
    ├── S_tortugarum.l.cairo.svg<br>
    ├── S_tortugarum.l.png<br>
    └── us.png<br>

7 directories, 177 files
```