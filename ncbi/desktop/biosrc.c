/*   biosrc.c
* ===========================================================================
*
*                            PUBLIC DOMAIN NOTICE
*            National Center for Biotechnology Information (NCBI)
*
*  This software/database is a "United States Government Work" under the
*  terms of the United States Copyright Act.  It was written as part of
*  the author's official duties as a United States Government employee and
*  thus cannot be copyrighted.  This software/database is freely available
*  to the public for use. The National Library of Medicine and the U.S.
*  Government do not place any restriction on its use or reproduction.
*  We would, however, appreciate having the NCBI and the author cited in
*  any work or product based on this material
*
*  Although all reasonable efforts have been taken to ensure the accuracy
*  and reliability of the software and data, the NLM and the U.S.
*  Government do not and cannot warrant the performance or results that
*  may be obtained by using this software or data. The NLM and the U.S.
*  Government disclaim all warranties, express or implied, including
*  warranties of performance, merchantability or fitness for any particular
*  purpose.
*
* ===========================================================================
*
* File Name:  biosrc.c
*
* Author:  Jonathan Kans
*
* Version Creation Date:   1/22/95
*
* $Revision: 6.93 $
*
* File Description: 
*
* Modifications:  
* --------------------------------------------------------------------------
* Date     Name        Description of modification
* -------  ----------  -----------------------------------------------------
*
*
* ==========================================================================
*/

#include <biosrc.h>
#include <document.h>
#include <gather.h>
#include <subutil.h>
#include <explore.h>
#define NLM_GENERATED_CODE_PROTO
#include <objmacro.h>
#include <macroapi.h>


extern EnumFieldAssoc  biosource_genome_simple_alist [];
ENUM_ALIST(biosource_genome_simple_alist)
  {" ",                    0},
  {"Genomic",              1},
  {"Chloroplast",          2},
  {"Kinetoplast",          4},
  {"Mitochondrion",        5},
  {"Plastid",              6},
  {"Macronuclear",         7},
  {"Extrachromosomal",     8},
  {"Plasmid",              9},
    /*
  {"Transposon",          10},
  {"Insertion Sequence",  11},
    */
  {"Cyanelle",            12},
  {"Proviral",            13},
  {"Virion",              14},
  {"Nucleomorph",         15},
  {"Apicoplast",          16},
  {"Leucoplast",          17},
  {"Proplastid",          18},
  {"Endogenous-virus",    19},
  {"Hydrogenosome",       20},
  {"Chromosome",          21},
  {"Chromatophore",       22},
END_ENUM_ALIST


extern EnumFieldAssoc  biosource_origin_alist [];
ENUM_ALIST(biosource_origin_alist)
  {" ",               0},
  {"Natural",         1},
  {"Natural Mutant",  2},
  {"Mutant",          3},
  {"Artificial",      4},
  {"Synthetic",       5},
  {"Other",         255},
END_ENUM_ALIST

Int2     numGeneticCodes = 0;
Int2     gcIdToIndex [NUM_GENETIC_CODES];
Uint1    gcIndexToId [NUM_GENETIC_CODES];
CharPtr  gcNames [NUM_GENETIC_CODES];

static CharPtr  orgTxtPtr = NULL;
static CharPtr  PNTR orgStrIdx = NULL;
static Int2     orgNum = 0;

#define ORGANISM_PAGE         0
#define MODIFIERS_PAGE        1
#define MISCELLANEOUS_PAGE    2
#define COMMON_PAGE           3
#define LOCATION_PAGE         4

#define NUM_PAGES  8

typedef struct genbiopage {
  DIALOG_MESSAGE_BLOCK
  TexT            taxName;
  Handle          commonName;
  Boolean         typedSciName;
  Boolean         typedComName;
  Int2            selectedOrg;
  Int2            clickedOrg;
  DoC             orglist;
  Int2            nuclGC;
  Int2            mitoGC;
  Int4            taxID;
  DialoG          genome;
  PopuP           origin;
  ButtoN          is_focus;
  PopuP           simplecode;
  PopuP           gcode;
  PopuP           mgcode;
  TexT            lineage;
  TexT            gbDiv;
  DialoG          db;
  DialoG          syn;
  DialoG          mod;
  GrouP           orgGrp [5];
  GrouP           modGrp [5];
  GrouP           miscGrp [5];
  PrompT          gbacr;
  PrompT          gbana;
  PrompT          gbsyn;

  DialoG          subsrc_val_dlg;
  DialoG          orgmod_val_dlg;

  CharPtr         origTaxName;
  Boolean         stripOldName;
//  EnumFieldAssoc  PNTR genomeAlist;
  Uint1           orgname_choice;
  Pointer         orgname_data;
  EnumFieldAssocPtr orgmod_alists [2];
  EnumFieldAssocPtr subsource_alists [2];
} GenBioPage, PNTR GenBioPagePtr;

typedef struct genbioform {
  FEATURE_FORM_BLOCK
  SeqEntryPtr   sep;
  GrouP         pages [NUM_PAGES];
  DialoG        foldertabs;
  Int2          currentPage;

  LookupTaxonomyProc  lookupTaxonomy;
} GenBioForm, PNTR GenBioFormPtr;

#ifndef WIN16
static CharPtr taxlistMemStrs [] = {
  "15\n",
  "Acanthamoeba castellanii\t\t1\t4\tINV\t5755\n",
  "Acanthoscurria gomesiana\t\t1\t5\tINV\t115339\n",
  "Acetabularia acetabulum\t\t6\t1\tPLN\t35845\n",
  "Acipenser sinensis\tChinese sturgeon\t1\t2\tVRT\t61970\n",
  "Acipenser transmontanus\twhite sturgeon\t1\t2\tVRT\t7904\n",
  "Acorus americanus\t\t1\t1\tPLN\t263995\n",
  "Acropora millepora\t\t1\t4\tINV\t45264\n",
  "Acropora palmata\t\t1\t4\tINV\t6131\n",
  "Acyrthosiphon pisum\tpea aphid\t1\t5\tINV\t7029\n",
  "Adiantum capillus-veneris\t\t1\t1\tPLN\t13818\n",
  "Aedes aegypti\tyellow fever mosquito\t1\t5\tINV\t7159\n",
  "Aegilops speltoides\t\t1\t1\tPLN\t4573\n",
  "Aegilops tauschii\t\t1\t1\tPLN\t37682\n",
  "Agrostis capillaris\t\t1\t1\tPLN\t204232\n",
  "Agrostis stolonifera\t\t1\t1\tPLN\t63632\n",
  "Ajellomyces capsulatus\t\t1\t4\tPLN\t5037\n",
  "Ajellomyces capsulatus NAm1\t\t1\t4\tPLN\t339724\n",
  "Alexandrium tamarense\t\t1\t4\tPLN\t2926\n",
  "Alligator mississippiensis\tAmerican alligator\t1\t2\tVRT\t8496\n",
  "Allium cepa\tonion\t1\t1\tPLN\t4679\n",
  "Allomyces macrogynus\t\t1\t4\tPLN\t28583\n",
  "Alternaria brassicicola\t\t1\t4\tPLN\t29001\n",
  "Amblyomma americanum\tlone star tick\t1\t5\tINV\t6943\n",
  "Amblyomma variegatum\t\t1\t5\tINV\t34610\n",
  "Amborella trichopoda\t\t1\t1\tPLN\t13333\n",
  "Ambystoma mexicanum\taxolotl\t1\t2\tVRT\t8296\n",
  "Ambystoma ordinarium\tPuerto Hondo stream salamander\t1\t2\tVRT\t288796\n",
  "Ambystoma tigrinum tigrinum\tEastern tiger salamander\t1\t2\tVRT\t43116\n",
  "Amoebidium parasiticum\t\t1\t4\tINV\t4881\n",
  "Amorphotheca resinae\tcreosote fungus\t1\t4\tPLN\t5101\n",
  "Amphidinium carterae\t\t1\t4\tPLN\t2961\n",
  "Ananas comosus\tpineapple\t1\t1\tPLN\t4615\n",
  "Anas platyrhynchos\t\t1\t2\tVRT\t8839\n",
  "Ancylostoma caninum\tdog hookworm\t1\t5\tINV\t29170\n",
  "Ancylostoma ceylanicum\t\t1\t5\tINV\t53326\n",
  "Anolis carolinensis\tgreen anole\t1\t2\tVRT\t28377\n",
  "Anolis sagrei\tbrown anole\t1\t2\tVRT\t38937\n",
  "Anopheles albimanus\t\t1\t5\tINV\t7167\n",
  "Anopheles funestus\tAfrican malaria mosquito\t1\t5\tINV\t62324\n",
  "Anopheles gambiae\tAfrican malaria mosquito\t1\t5\tINV\t7165\n",
  "Anopheles gambiae str. PEST\t\t1\t5\tINV\t180454\n",
  "Antheraea mylitta\t\t1\t5\tINV\t34739\n",
  "Antirrhinum majus\tsnapdragon\t1\t1\tPLN\t4151\n",
  "Antonospora locustae\t\t1\t1\tINV\t278021\n",
  "Antrodia cinnamomea\t\t1\t4\tPLN\t279009\n",
  "Aphanomyces cochlioides\t\t1\t1\tPLN\t112091\n",
  "Aphis gossypii\tcotton aphid\t1\t5\tINV\t80765\n",
  "Apis mellifera\thoney bee\t1\t5\tINV\t7460\n",
  "Aplysia californica\tCalifornia sea hare\t1\t5\tINV\t6500\n",
  "Aquilegia formosa\t\t1\t1\tPLN\t223430\n",
  "Arabidopsis lyrata subsp. petraea\t\t1\t1\tPLN\t59691\n",
  "Arabidopsis thaliana\tthale cress\t1\t1\tPLN\t3702\n",
  "Arachis batizocoi\t\t1\t1\tPLN\t108210\n",
  "Arachis duranensis\t\t1\t1\tPLN\t130453\n",
  "Arachis hypogaea\tpeanut\t1\t1\tPLN\t3818\n",
  "Arachis stenosperma\t\t1\t1\tPLN\t217475\n",
  "Argas monolakensis\t\t1\t5\tINV\t34602\n",
  "Argopecten irradians\t\t1\t5\tINV\t31199\n",
  "Ascaris suum\tpig roundworm\t1\t5\tINV\t6253\n",
  "Ascosphaera apis USDA-ARSEF 7405\t\t1\t4\tPLN\t392613\n",
  "Ashbya gossypii ATCC 10895\t\t1\t3\tPLN\t284811\n",
  "Asparagus officinalis\tgarden asparagus\t1\t1\tPLN\t4686\n",
  "Aspergillus clavatus NRRL 1\t\t1\t4\tPLN\t344612\n",
  "Aspergillus flavus\t\t1\t4\tPLN\t5059\n",
  "Aspergillus flavus NRRL3357\t\t1\t4\tPLN\t332952\n",
  "Aspergillus fumigatus Af293\t\t1\t4\tPLN\t330879\n",
  "Aspergillus nidulans FGSC A4\t\t1\t4\tPLN\t227321\n",
  "Aspergillus niger\t\t1\t4\tPLN\t5061\n",
  "Aspergillus niger CBS 513.88\t\t1\t4\tPLN\t425011\n",
  "Aspergillus oryzae\t\t1\t4\tPLN\t5062\n",
  "Aspergillus terreus NIH2624\t\t1\t4\tPLN\t341663\n",
  "Astatotilapia burtoni\t\t1\t2\tVRT\t8153\n",
  "Aureobasidium pullulans\t\t1\t4\tPLN\t5580\n",
  "Avena sativa\toat\t1\t1\tPLN\t4498\n",
  "Babesia bovis\t\t1\t4\tINV\t5865\n",
  "Bacillus cereus\t\t11\t0\tBCT\t1396\n",
  "Bacillus clausii\t\t11\t0\tBCT\t79880\n",
  "Bacillus licheniformis\t\t11\t0\tBCT\t1402\n",
  "Bacillus subtilis\t\t11\t0\tBCT\t1423\n",
  "Bemisia tabaci\t\t1\t5\tINV\t7038\n",
  "Beta vulgaris\t\t1\t1\tPLN\t161934\n",
  "Betula pendula\tEuropean white birch\t1\t1\tPLN\t3505\n",
  "Bicyclus anynana\tsquinting bush brown\t1\t5\tINV\t110368\n",
  "Bigelowiella natans\t\t1\t1\tINV\t227086\n",
  "Biomphalaria glabrata\t\t1\t5\tINV\t6526\n",
  "Blastocladiella emersonii\t\t1\t4\tPLN\t4808\n",
  "Blastocystis hominis\t\t1\t1\tPLN\t12968\n",
  "Blumeria graminis f. sp. hordei\t\t1\t4\tPLN\t62688\n",
  "Boechera stricta\t\t1\t1\tPLN\t72658\n",
  "Bombyx mori\tdomestic silkworm\t1\t5\tINV\t7091\n",
  "Borrelia burgdorferi\tLyme disease spirochete\t11\t0\tBCT\t139\n",
  "Bos indicus\t\t1\t2\tMAM\t9915\n",
  "Bos taurus\tcattle\t1\t2\tMAM\t9913\n",
  "Botryotinia fuckeliana\t\t1\t4\tPLN\t40559\n",
  "Botryotinia fuckeliana B05.10\t\t1\t4\tPLN\t332648\n",
  "Brachionus plicatilis\t\t1\t5\tINV\t10195\n",
  "Brachypodium distachyon\t\t1\t1\tPLN\t15368\n",
  "Bradyrhizobium japonicum\t\t11\t0\tBCT\t375\n",
  "Branchiostoma floridae\tFlorida lancelet\t1\t5\tINV\t7739\n",
  "Brassica carinata\t\t1\t1\tPLN\t52824\n",
  "Brassica napus\trape\t1\t1\tPLN\t3708\n",
  "Brassica oleracea\t\t1\t1\tPLN\t3712\n",
  "Brassica oleracea var. alboglabra\tChinese kale\t1\t1\tPLN\t3714\n",
  "Brassica oleracea var. italica\tasparagus broccoli\t1\t1\tPLN\t36774\n",
  "Brassica rapa\t\t1\t1\tPLN\t3711\n",
  "Brassica rapa subsp. pekinensis\t\t1\t1\tPLN\t51351\n",
  "Brucella abortus\t\t11\t0\tBCT\t235\n",
  "Brugia malayi\t\t1\t5\tINV\t6279\n",
  "Bruguiera gymnorhiza\t\t1\t1\tPLN\t39984\n",
  "Bubalus bubalis\twater buffalo\t1\t2\tMAM\t89462\n",
  "Burkholderia pseudomallei 112\t\t11\t0\tBCT\t441154\n",
  "Burkholderia pseudomallei 14\t\t11\t0\tBCT\t441160\n",
  "Burkholderia pseudomallei 381\t\t11\t0\tBCT\t441157\n",
  "Burkholderia pseudomallei 7894\t\t11\t0\tBCT\t441156\n",
  "Burkholderia pseudomallei 9\t\t11\t0\tBCT\t441158\n",
  "Burkholderia pseudomallei 91\t\t11\t0\tBCT\t441159\n",
  "Burkholderia pseudomallei B7210\t\t11\t0\tBCT\t441155\n",
  "Burkholderia pseudomallei DM98\t\t11\t0\tBCT\t441161\n",
  "Bursaphelenchus mucronatus\t\t1\t5\tINV\t6325\n",
  "Bursaphelenchus xylophilus\t\t1\t5\tINV\t6326\n",
  "Caenorhabditis brenneri\t\t1\t5\tINV\t135651\n",
  "Caenorhabditis briggsae\t\t1\t5\tINV\t6238\n",
  "Caenorhabditis briggsae AF16\t\t1\t5\tINV\t473542\n",
  "Caenorhabditis elegans\t\t1\t5\tINV\t6239\n",
  "Caenorhabditis remanei\t\t1\t5\tINV\t31234\n",
  "Calanus finmarchicus\t\t1\t5\tINV\t6837\n",
  "Callinectes sapidus\tblue crab\t1\t5\tINV\t6763\n",
  "Callithrix jacchus\twhite-tufted-ear marmoset\t1\t2\tPRI\t9483\n",
  "Callorhinchus milii\telephantfish\t1\t2\tVRT\t7868\n",
  "Camellia sinensis\t\t1\t1\tPLN\t4442\n",
  "Candida albicans\t\t12\t4\tPLN\t5476\n",
  "Candida albicans SC5314\t\t12\t4\tPLN\t237561\n",
  "Candida glabrata\t\t1\t3\tPLN\t5478\n",
  "Candida glabrata CBS 138\t\t1\t3\tPLN\t284593\n",
  "Candida parapsilosis\t\t12\t4\tPLN\t5480\n",
  "Candida tropicalis\t\t12\t3\tPLN\t5482\n",
  "candidate division TM7 single-cell isolate TM7a\t\t11\t0\tBCT\t447454\n",
  "Canis latrans\tcoyote\t1\t2\tMAM\t9614\n",
  "Canis lupus\tgray wolf\t1\t2\tMAM\t9612\n",
  "Canis lupus familiaris\tdog\t1\t2\tMAM\t9615\n",
  "Cannabis sativa\themp\t1\t1\tPLN\t3483\n",
  "Capra hircus\tgoat\t1\t2\tMAM\t9925\n",
  "Capsaspora owczarzaki\t\t1\t1\tINV\t192875\n",
  "Capsicum annuum\t\t1\t1\tPLN\t4072\n",
  "Carcinus maenas\tgreen crab\t1\t5\tINV\t6759\n",
  "Carica papaya\tpapaya\t1\t1\tPLN\t3649\n",
  "Carthamus tinctorius\tsafflower\t1\t1\tPLN\t4222\n",
  "Catharanthus roseus\tMadagascar periwinkle\t1\t1\tPLN\t4058\n",
  "Cavia porcellus\tdomestic guinea pig\t1\t2\tROD\t10141\n",
  "Celuca pugilator\tAtlantic sand fiddler crab\t1\t5\tINV\t6772\n",
  "Cenchrus ciliaris\t\t1\t1\tPLN\t35872\n",
  "Centaurea maculosa\t\t1\t1\tPLN\t215693\n",
  "Centaurea solstitialis\t\t1\t1\tPLN\t347529\n",
  "Ceratodon purpureus\t\t1\t1\tPLN\t3225\n",
  "Ceratopteris richardii\t\t1\t1\tPLN\t49495\n",
  "Cercomonas longicauda\t\t1\t1\tINV\t100933\n",
  "Chaetomium cupreum\t\t1\t4\tPLN\t155874\n",
  "Chaetomium globosum CBS 148.51\t\t1\t4\tPLN\t306901\n",
  "Chamaecyparis obtusa\t\t1\t1\tPLN\t13415\n",
  "Chironomus tentans\t\t1\t5\tINV\t7153\n",
  "Chlamydia trachomatis\t\t11\t0\tBCT\t813\n",
  "Chlamydomonas incerta\t\t1\t1\tPLN\t51695\n",
  "Chlamydomonas reinhardtii\t\t1\t1\tPLN\t3055\n",
  "Chlamys farreri\t\t1\t5\tINV\t202578\n",
  "Chlorocebus aethiops\tAfrican green monkey\t1\t2\tPRI\t9534\n",
  "Chondrus crispus\tcarragheen\t1\t4\tPLN\t2769\n",
  "Chrysemys picta\t\t1\t2\tVRT\t8479\n",
  "Cicer arietinum\tchickpea\t1\t1\tPLN\t3827\n",
  "Cichorium endivia\t\t1\t1\tPLN\t114280\n",
  "Cichorium intybus\tchicory\t1\t1\tPLN\t13427\n",
  "Ciona intestinalis\t\t1\t13\tINV\t7719\n",
  "Ciona savignyi\t\t1\t13\tINV\t51511\n",
  "Citrus aurantium\t\t1\t1\tPLN\t43166\n",
  "Citrus clementina\t\t1\t1\tPLN\t85681\n",
  "Citrus reshni\t\t1\t1\tPLN\t171252\n",
  "Citrus reticulata\t\t1\t1\tPLN\t85571\n",
  "Citrus sinensis\t\t1\t1\tPLN\t2711\n",
  "Citrus unshiu\t\t1\t1\tPLN\t55188\n",
  "Cleome hassleriana\t\t1\t1\tPLN\t28532\n",
  "Clonorchis sinensis\t\t1\t9\tINV\t79923\n",
  "Closterium peracerosum-strigosum-littorale complex\t\t1\t1\tPLN\t34146\n",
  "Coccidioides immitis H538.4\t\t1\t4\tPLN\t396776\n",
  "Coccidioides immitis RMSCC 3703\t\t1\t4\tPLN\t454286\n",
  "Coccidioides immitis RS\t\t1\t4\tPLN\t246410\n",
  "Coccidioides posadasii\t\t1\t4\tPLN\t199306\n",
  "Coccidioides posadasii str. Silveira\t\t1\t4\tPLN\t443226\n",
  "Coffea arabica\tcoffee\t1\t1\tPLN\t13443\n",
  "Coffea canephora\t\t1\t1\tPLN\t49390\n",
  "Coprinopsis cinerea\t\t1\t4\tPLN\t5346\n",
  "Cordyceps bassiana\t\t1\t4\tPLN\t176275\n",
  "Coregonus clupeaformis\tlake whitefish\t1\t2\tVRT\t59861\n",
  "Corynascus heterothallicus\t\t1\t4\tPLN\t78579\n",
  "Corynebacterium glutamicum\t\t11\t0\tBCT\t1718\n",
  "Crassostrea gigas\tPacific oyster\t1\t5\tINV\t29159\n",
  "Crassostrea virginica\teastern oyster\t1\t5\tINV\t6565\n",
  "Cricetulus griseus\tChinese hamster\t1\t2\tROD\t10029\n",
  "Crocus sativus\t\t1\t1\tPLN\t82528\n",
  "Cryphonectria parasitica\t\t1\t4\tPLN\t5116\n",
  "Cryptococcus laurentii\t\t1\t4\tPLN\t5418\n",
  "Cryptococcus neoformans var. neoformans\t\t1\t4\tPLN\t40410\n",
  "Cryptococcus neoformans var. neoformans B-3501A\t\t1\t4\tPLN\t283643\n",
  "Cryptococcus neoformans var. neoformans JEC21\t\t1\t4\tPLN\t214684\n",
  "Cryptomeria japonica\tJapanese cedar\t1\t1\tPLN\t3369\n",
  "Cryptosporidium hominis TU502\t\t1\t4\tINV\t353151\n",
  "Cryptosporidium parvum\t\t1\t4\tINV\t5807\n",
  "Cryptosporidium parvum Iowa II\t\t1\t4\tINV\t353152\n",
  "Ctenocephalides felis\tcat flea\t1\t5\tINV\t7515\n",
  "Cucumis melo\tmuskmelon\t1\t1\tPLN\t3656\n",
  "Cucumis melo subsp. agrestis\t\t1\t1\tPLN\t217619\n",
  "Cucumis melo subsp. melo\t\t1\t1\tPLN\t412675\n",
  "Cucumis sativus\tcucumber\t1\t1\tPLN\t3659\n",
  "Culex pipiens quinquefasciatus\tsouthern house mosquito\t1\t5\tINV\t7176\n",
  "Culicoides sonorensis\t\t1\t5\tINV\t179676\n",
  "Cunninghamella elegans\t\t1\t4\tPLN\t4853\n",
  "Curcuma longa\t\t1\t1\tPLN\t136217\n",
  "Cyamopsis tetragonoloba\tguar\t1\t1\tPLN\t3832\n",
  "Cyanidioschyzon merolae strain 10D\t\t1\t1\tPLN\t280699\n",
  "Cyanophora paradoxa\t\t1\t1\tPLN\t2762\n",
  "Cycas rumphii\t\t1\t1\tPLN\t58031\n",
  "Cynodon dactylon\tBermuda grass\t1\t1\tPLN\t28909\n",
  "Cyprinus carpio\tcommon carp\t1\t2\tVRT\t7962\n",
  "Danio rerio\tzebrafish\t1\t2\tVRT\t7955\n",
  "Daphnia magna\t\t1\t5\tINV\t35525\n",
  "Dasypus novemcinctus\tnine-banded armadillo\t1\t2\tMAM\t9361\n",
  "Debaryomyces hansenii CBS767\t\t1\t3\tPLN\t284592\n",
  "Debaryomyces hansenii var. hansenii\t\t1\t3\tPLN\t58641\n",
  "Deinagkistrodon acutus\t\t1\t2\tVRT\t36307\n",
  "Dekkera bruxellensis\t\t1\t4\tPLN\t5007\n",
  "Diabrotica virgifera virgifera\twestern corn rootworm\t1\t5\tINV\t50390\n",
  "Diaphorina citri\tAsian citrus psyllid\t1\t5\tINV\t121845\n",
  "Diaprepes abbreviatus\tDiaprepes root weevil\t1\t5\tINV\t13040\n",
  "Dicentrarchus labrax\tEuropean sea bass\t1\t2\tVRT\t13489\n",
  "Dictyocaulus viviparus\tbovine lungworm\t1\t5\tINV\t29172\n",
  "Dictyostelium discoideum\t\t1\t1\tINV\t44689\n",
  "Dictyostelium discoideum AX4\t\t1\t1\tINV\t352472\n",
  "Diplonema papillatum\t\t1\t4\tINV\t91374\n",
  "Dirofilaria immitis\tdog heartworm nematode\t1\t5\tINV\t6287\n",
  "Drosophila ananassae\t\t1\t5\tINV\t7217\n",
  "Drosophila auraria\t\t1\t5\tINV\t47315\n",
  "Drosophila erecta\t\t1\t5\tINV\t7220\n",
  "Drosophila grimshawi\t\t1\t5\tINV\t7222\n",
  "Drosophila melanogaster\tfruit fly\t1\t5\tINV\t7227\n",
  "Drosophila mojavensis\t\t1\t5\tINV\t7230\n",
  "Drosophila persimilis\t\t1\t5\tINV\t7234\n",
  "Drosophila pseudoobscura\t\t1\t5\tINV\t7237\n",
  "Drosophila sechellia\t\t1\t5\tINV\t7238\n",
  "Drosophila simulans\t\t1\t5\tINV\t7240\n",
  "Drosophila virilis\t\t1\t5\tINV\t7244\n",
  "Drosophila willistoni\t\t1\t5\tINV\t7260\n",
  "Drosophila yakuba\t\t1\t5\tINV\t7245\n",
  "Dugesia japonica\t\t1\t9\tINV\t6161\n",
  "Dugesia ryukyuensis\t\t1\t9\tINV\t79738\n",
  "Dunaliella salina\t\t1\t1\tPLN\t3046\n",
  "Echinococcus granulosus\t\t1\t9\tINV\t6210\n",
  "Echinops telfairi\tsmall Madagascar hedgehog\t1\t2\tMAM\t9371\n",
  "Eimeria tenella\t\t1\t4\tINV\t5802\n",
  "Eisenia fetida\tcommon brandling worm\t1\t5\tINV\t6396\n",
  "Elaeis guineensis\tAfrican oil palm\t1\t1\tPLN\t51953\n",
  "Elaeis oleifera\t\t1\t1\tPLN\t80265\n",
  "Emericella nidulans\t\t1\t4\tPLN\t162425\n",
  "Emiliania huxleyi\t\t1\t4\tPLN\t2903\n",
  "Endoriftia persephone 'Hot96_1+Hot96_2'\t\t11\t0\tBCT\t394104\n",
  "Entamoeba dispar\t\t1\t1\tINV\t46681\n",
  "Entamoeba dispar SAW760\t\t1\t1\tINV\t370354\n",
  "Entamoeba histolytica\t\t1\t1\tINV\t5759\n",
  "Entamoeba histolytica HM-1:IMSS\t\t1\t1\tINV\t294381\n",
  "Entamoeba invadens\t\t1\t1\tINV\t33085\n",
  "Entamoeba invadens IP1\t\t1\t1\tINV\t370355\n",
  "Entamoeba moshkovskii FIC\t\t1\t1\tINV\t434306\n",
  "Entamoeba terrapinae\t\t1\t1\tINV\t110771\n",
  "Epiphyas postvittana\t\t1\t5\tINV\t65032\n",
  "Eptatretus burgeri\tinshore hagfish\t1\t2\tVRT\t7764\n",
  "Equus caballus\thorse\t1\t2\tMAM\t9796\n",
  "Eragrostis tef\t\t1\t1\tPLN\t110835\n",
  "Erinaceus europaeus\twestern European hedgehog\t1\t2\tMAM\t9365\n",
  "Escherichia coli\t\t11\t0\tBCT\t562\n",
  "Eschscholzia californica\tCalifornia poppy\t1\t1\tPLN\t3467\n",
  "Esox lucius\tnorthern pike\t1\t2\tVRT\t8010\n",
  "Eubalaena glacialis\tNorth Atlantic right whale\t1\t2\tMAM\t27606\n",
  "Eucalyptus camaldulensis\tMurray red gum\t1\t1\tPLN\t34316\n",
  "Eucalyptus globulus subsp. bicostata\t\t1\t1\tPLN\t71272\n",
  "Eucalyptus grandis\t\t1\t1\tPLN\t71139\n",
  "Eucalyptus gunnii\tcider tree\t1\t1\tPLN\t3933\n",
  "Euglena gracilis\t\t1\t4\tPLN\t3039\n",
  "Euglena longa\t\t1\t4\tPLN\t3037\n",
  "Euphorbia esula\tleafy spurge\t1\t1\tPLN\t3993\n",
  "Euphorbia tirucalli\t\t1\t1\tPLN\t142860\n",
  "Euprymna scolopes\t\t1\t5\tINV\t6613\n",
  "Felis catus\tdomestic cat\t1\t2\tMAM\t9685\n",
  "Fenneropenaeus chinensis\t\t1\t5\tINV\t139456\n",
  "Festuca arundinacea\t\t1\t1\tPLN\t4606\n",
  "Folsomia candida\t\t1\t5\tINV\t158441\n",
  "fossil metagenome\t\t11\t2\tENV\t444079\n",
  "Fragaria vesca\t\t1\t1\tPLN\t57918\n",
  "Fragaria vesca subsp. vesca\t\t1\t1\tPLN\t101020\n",
  "Fragilariopsis cylindrus\t\t1\t1\tPLN\t186039\n",
  "Fundulus heteroclitus\tkillifish\t1\t2\tVRT\t8078\n",
  "Fusarium oxysporum f. sp. cucumerinum\t\t1\t4\tPLN\t5508\n",
  "Fusarium oxysporum f. sp. melonis\t\t1\t4\tPLN\t61369\n",
  "Fusarium solani\t\t1\t4\tPLN\t169388\n",
  "Fusarium sporotrichioides\t\t1\t4\tPLN\t5514\n",
  "Fusarium virguliforme\t\t1\t4\tPLN\t232082\n",
  "Gadus morhua\tAtlantic cod\t1\t2\tVRT\t8049\n",
  "Gallus gallus\tchicken\t1\t2\tVRT\t9031\n",
  "Gammarus pulex\t\t1\t5\tINV\t52641\n",
  "Gasterosteus aculeatus\tthree-spined stickleback\t1\t2\tVRT\t69293\n",
  "GB virus C\t\t1\t0\tVRL\t54290\n",
  "Gekko japonicus\t\t1\t2\tVRT\t146911\n",
  "Gene trapping vector VICTR75\t\t11\t0\tSYN\t447634\n",
  "Gene trapping vector VICTR76\t\t11\t0\tSYN\t447635\n",
  "Geomyces pannorum\t\t1\t4\tPLN\t79858\n",
  "Gerbera hybrid cv. 'Terra Regina'\t\t1\t1\tPLN\t226891\n",
  "Giardia intestinalis\t\t1\t0\tINV\t5741\n",
  "Giardia lamblia ATCC 50803\t\t1\t0\tINV\t184922\n",
  "Gibberella moniliformis\t\t1\t4\tPLN\t117187\n",
  "Gibberella zeae\t\t1\t4\tPLN\t5518\n",
  "Gibberella zeae PH-1\t\t1\t4\tPLN\t229533\n",
  "Gillichthys mirabilis\tlong-jawed mudsucker\t1\t2\tVRT\t8222\n",
  "Ginkgo biloba\tmaidenhair tree\t1\t1\tPLN\t3311\n",
  "Glaucocystis nostochinearum\t\t1\t1\tPLN\t38271\n",
  "Globodera pallida\t\t1\t5\tINV\t36090\n",
  "Globodera rostochiensis\t\t1\t5\tINV\t31243\n",
  "Gloeophyllum trabeum\t\t1\t4\tPLN\t104355\n",
  "Glomus intraradices\t\t1\t4\tPLN\t4876\n",
  "Glossina morsitans morsitans\t\t1\t5\tINV\t37546\n",
  "Glycine max\tsoybean\t1\t1\tPLN\t3847\n",
  "Glycine soja\t\t1\t1\tPLN\t3848\n",
  "Glycyphagus domesticus\t\t1\t5\tINV\t105145\n",
  "Gnetum gnemon\t\t1\t1\tPLN\t3382\n",
  "Gobiocypris rarus\t\t1\t2\tVRT\t143606\n",
  "Gorilla gorilla\tgorilla\t1\t2\tPRI\t9593\n",
  "Gossypium arboreum\t\t1\t1\tPLN\t29729\n",
  "Gossypium exiguum\t\t1\t1\tPLN\t47626\n",
  "Gossypium herbaceum\t\t1\t1\tPLN\t34274\n",
  "Gossypium hirsutum\tupland cotton\t1\t1\tPLN\t3635\n",
  "Gossypium raimondii\t\t1\t1\tPLN\t29730\n",
  "Gracilaria changii\t\t1\t4\tPLN\t172969\n",
  "Graphocephala atropunctata\t\t1\t5\tINV\t36148\n",
  "Grosmannia clavigera\t\t1\t4\tPLN\t226899\n",
  "Gryllus bimaculatus\ttwo-spotted cricket\t1\t5\tINV\t6999\n",
  "Guillardia theta\t\t1\t1\tPLN\t55529\n",
  "Haemonchus contortus\t\t1\t5\tINV\t6289\n",
  "Haemophilus influenzae\t\t11\t0\tBCT\t727\n",
  "Halocynthia roretzi\t\t1\t13\tINV\t7729\n",
  "Hartmannella vermiformis\t\t1\t1\tINV\t5778\n",
  "Hebeloma cylindrosporum\t\t1\t4\tPLN\t76867\n",
  "Hedyotis centranthoides\t\t1\t1\tPLN\t219666\n",
  "Hedyotis terminalis\t\t1\t1\tPLN\t219667\n",
  "Helianthus annuus\tcommon sunflower\t1\t1\tPLN\t4232\n",
  "Helianthus argophyllus\t\t1\t1\tPLN\t73275\n",
  "Helianthus ciliaris\t\t1\t1\tPLN\t73280\n",
  "Helianthus exilis\t\t1\t1\tPLN\t400408\n",
  "Helianthus paradoxus\t\t1\t1\tPLN\t73304\n",
  "Helianthus petiolaris\t\t1\t1\tPLN\t4234\n",
  "Helianthus tuberosus\t\t1\t1\tPLN\t4233\n",
  "Helicobacter pylori\t\t11\t0\tBCT\t210\n",
  "Heliconius erato\tcrimson-patched longwing\t1\t5\tINV\t33431\n",
  "Heliconius melpomene\t\t1\t5\tINV\t34740\n",
  "Hemicentrotus pulcherrimus\t\t1\t9\tINV\t7650\n",
  "Hepatitis A virus\t\t1\t0\tVRL\t12092\n",
  "Hepatitis B virus\t\t1\t0\tVRL\t10407\n",
  "Hepatitis C virus\t\t1\t0\tVRL\t11103\n",
  "Hepatitis C virus subtype 1a\t\t1\t0\tVRL\t31646\n",
  "Hepatitis C virus subtype 1b\t\t1\t0\tVRL\t31647\n",
  "Heterocapsa triquetra\t\t1\t4\tPLN\t66468\n",
  "Heterodera glycines\t\t1\t5\tINV\t51029\n",
  "Heterodera schachtii\t\t1\t5\tINV\t97005\n",
  "Heterorhabditis bacteriophora\t\t1\t5\tINV\t37862\n",
  "Hevea brasiliensis\t\t1\t1\tPLN\t3981\n",
  "Hippoglossus hippoglossus\tAtlantic halibut\t1\t2\tVRT\t8267\n",
  "Histiona aroides\t\t1\t1\tINV\t392300\n",
  "Holothuria glaberrima\t\t1\t9\tINV\t31192\n",
  "Homalodisca coagulata\tglassy-winged sharpshooter\t1\t5\tINV\t197043\n",
  "Homarus americanus\tAmerican lobster\t1\t5\tINV\t6706\n",
  "Homo sapiens\thuman\t1\t2\tPRI\t9606\n",
  "Hordeum vulgare\t\t1\t1\tPLN\t4513\n",
  "Hordeum vulgare subsp. spontaneum\t\t1\t1\tPLN\t77009\n",
  "Hordeum vulgare subsp. vulgare\tdomesticated barley\t1\t1\tPLN\t112509\n",
  "Human herpesvirus 5\tHuman cytomegalovirus\t1\t0\tVRL\t10359\n",
  "Human immunodeficiency virus 1\t\t1\t0\tVRL\t11676\n",
  "Human immunodeficiency virus 2\t\t1\t0\tVRL\t11709\n",
  "Human metapneumovirus\t\t1\t0\tVRL\t162145\n",
  "Human T-lymphotropic virus 1\t\t1\t0\tVRL\t11908\n",
  "Humulus lupulus\tEuropean hop\t1\t1\tPLN\t3486\n",
  "Hydra magnipapillata\t\t1\t4\tINV\t6085\n",
  "Hydra vulgaris\t\t1\t4\tINV\t6087\n",
  "Hydractinia echinata\t\t1\t4\tINV\t35630\n",
  "Hyperamoeba dachnaya\t\t1\t1\tINV\t181200\n",
  "Hypocrea jecorina\t\t1\t4\tPLN\t51453\n",
  "Hypocrea lixii\t\t1\t4\tPLN\t5544\n",
  "Hypocrea virens\t\t1\t4\tPLN\t29875\n",
  "Hypsibius dujardini\t\t1\t5\tINV\t232323\n",
  "Ichthyophthirius multifiliis\t\t6\t4\tINV\t5932\n",
  "Ictalurus furcatus\tblue catfish\t1\t2\tVRT\t66913\n",
  "Ictalurus punctatus\tchannel catfish\t1\t2\tVRT\t7998\n",
  "Idiosepius paradoxus\t\t1\t5\tINV\t294707\n",
  "Ipomoea batatas\tsweet potato\t1\t1\tPLN\t4120\n",
  "Ipomoea nil\tJapanese morning glory\t1\t1\tPLN\t35883\n",
  "Isochrysis galbana\t\t1\t1\tPLN\t37099\n",
  "Ixodes scapularis\tblack-legged tick\t1\t5\tINV\t6945\n",
  "Jakoba bahamiensis\t\t1\t1\tINV\t221721\n",
  "Jakoba libera\t\t1\t1\tINV\t143017\n",
  "Juglans regia\tEnglish walnut\t1\t1\tPLN\t51240\n",
  "Karenia brevis\t\t1\t4\tPLN\t156230\n",
  "Karlodinium micrum\t\t1\t4\tPLN\t342587\n",
  "Kazachstania exigua\t\t1\t3\tPLN\t34358\n",
  "Kazachstania unispora\t\t1\t3\tPLN\t27294\n",
  "Kluyveromyces lactis\t\t1\t4\tPLN\t28985\n",
  "Kluyveromyces lactis NRRL Y-1140\t\t1\t4\tPLN\t284590\n",
  "Kluyveromyces marxianus\t\t1\t3\tPLN\t4911\n",
  "Kluyveromyces thermotolerans\t\t1\t3\tPLN\t4916\n",
  "Laccaria bicolor\t\t1\t4\tPLN\t29883\n",
  "Lachancea kluyveri\t\t1\t3\tPLN\t4934\n",
  "Lactuca perennis\t\t1\t1\tPLN\t43195\n",
  "Lactuca saligna\t\t1\t1\tPLN\t75948\n",
  "Lactuca sativa\t\t1\t1\tPLN\t4236\n",
  "Lactuca serriola\t\t1\t1\tPLN\t75943\n",
  "Lactuca virosa\t\t1\t1\tPLN\t75947\n",
  "Lama pacos\talpaca\t1\t2\tMAM\t30538\n",
  "Laminaria digitata\t\t1\t1\tPLN\t80365\n",
  "Laupala kohalensis\t\t1\t5\tINV\t109027\n",
  "Lawsonia intracellularis\t\t11\t0\tBCT\t29546\n",
  "Leishmania braziliensis\t\t1\t4\tINV\t5660\n",
  "Leishmania donovani chagasi\t\t1\t4\tINV\t44271\n",
  "Leishmania infantum JPCM5\t\t1\t4\tINV\t435258\n",
  "Leishmania major\t\t1\t4\tINV\t5664\n",
  "Leishmania major strain Friedlin\t\t1\t4\tINV\t347515\n",
  "Lentinula edodes\tshiitake mushroom\t1\t4\tPLN\t5353\n",
  "Lepeophtheirus salmonis\tsalmon louse\t1\t5\tINV\t72036\n",
  "Leptinotarsa decemlineata\tColorado potato beetle\t1\t5\tINV\t7539\n",
  "Leucoraja erinacea\tlittle skate\t1\t2\tVRT\t7782\n",
  "Leucosporidium scottii\t\t1\t4\tPLN\t5278\n",
  "Limnanthes alba subsp. versicolor\t\t1\t1\tPLN\t377281\n",
  "Limonium bicolor\t\t1\t1\tPLN\t293754\n",
  "Lingulodinium polyedrum\t\t1\t4\tPLN\t160621\n",
  "Linum usitatissimum\tflax\t1\t1\tPLN\t4006\n",
  "Liriodendron tulipifera\t\t1\t1\tPLN\t3415\n",
  "Listeria monocytogenes\t\t11\t0\tBCT\t1639\n",
  "Listeria monocytogenes FSL F2-515\t\t11\t0\tBCT\t393116\n",
  "Listeria monocytogenes FSL J1-208\t\t11\t0\tBCT\t393119\n",
  "Lithognathus mormyrus\t\t1\t2\tVRT\t50591\n",
  "Litomosoides sigmodontis\t\t1\t5\tINV\t42156\n",
  "Litopenaeus vannamei\tPacific white shrimp\t1\t5\tINV\t6689\n",
  "Locusta migratoria\tmigratory locust\t1\t5\tINV\t7004\n",
  "Lodderomyces elongisporus NRRL YB-4239\t\t1\t3\tPLN\t379508\n",
  "Lolium multiflorum\tItalian ryegrass\t1\t1\tPLN\t4521\n",
  "Lolium temulentum\t\t1\t1\tPLN\t34176\n",
  "Lotus japonicus\t\t1\t1\tPLN\t34305\n",
  "Loxodonta africana\tAfrican savanna elephant\t1\t2\tMAM\t9785\n",
  "Lumbricus rubellus\thumus earthworm\t1\t5\tINV\t35632\n",
  "Lupinus albus\twhite lupine\t1\t1\tPLN\t3870\n",
  "Lutzomyia longipalpis\t\t1\t5\tINV\t7200\n",
  "Lycoris longituba\t\t1\t1\tPLN\t272140\n",
  "Lysiphlebus testaceipes\t\t1\t5\tINV\t77504\n",
  "Macaca fascicularis\tcrab-eating macaque\t1\t2\tPRI\t9541\n",
  "Macaca mulatta\trhesus monkey\t1\t2\tPRI\t9544\n",
  "Macaca nemestrina\tpig-tailed macaque\t1\t2\tPRI\t9545\n",
  "Maconellicoccus hirsutus\thibiscus mealybug\t1\t5\tINV\t177089\n",
  "Macropus eugenii\ttammar wallaby\t1\t2\tMAM\t9315\n",
  "Macrostomum lignano\t\t1\t9\tINV\t282301\n",
  "Magnaporthe grisea\t\t1\t4\tPLN\t148305\n",
  "Magnaporthe grisea 70-15\t\t1\t4\tPLN\t242507\n",
  "Magnetospirillum magnetotacticum MS-1\t\t11\t0\tBCT\t272627\n",
  "Malawimonas californiana\t\t1\t1\tINV\t221722\n",
  "Malawimonas jakobiformis\t\t1\t1\tINV\t136089\n",
  "Manduca sexta\ttobacco hornworm\t1\t5\tINV\t7130\n",
  "Manihot esculenta\tcassava\t1\t1\tPLN\t3983\n",
  "Marchantia polymorpha\tliverwort\t1\t1\tPLN\t3197\n",
  "marine metagenome\t\t11\t2\tENV\t408172\n",
  "Marsupenaeus japonicus\t\t1\t5\tINV\t27405\n",
  "Mastigamoeba balamuthi\t\t1\t1\tINV\t108607\n",
  "Measles virus\t\t1\t0\tVRL\t11234\n",
  "Medicago sativa\t\t1\t1\tPLN\t3879\n",
  "Medicago truncatula\tbarrel medic\t1\t1\tPLN\t3880\n",
  "Meleagris gallopavo\tturkey\t1\t2\tVRT\t9103\n",
  "Meloidogyne arenaria\t\t1\t5\tINV\t6304\n",
  "Meloidogyne chitwoodi\t\t1\t5\tINV\t59747\n",
  "Meloidogyne hapla\t\t1\t5\tINV\t6305\n",
  "Meloidogyne incognita\tsouthern root-knot nematode\t1\t5\tINV\t6306\n",
  "Meloidogyne javanica\troot-knot nematode\t1\t5\tINV\t6303\n",
  "Meloidogyne paranaensis\t\t1\t5\tINV\t189293\n",
  "Mesembryanthemum crystallinum\tcommon iceplant\t1\t1\tPLN\t3544\n",
  "Mesocricetus auratus\tgolden hamster\t1\t2\tROD\t10036\n",
  "Mesostigma viride\t\t1\t1\tPLN\t41882\n",
  "metagenome sequence\t\t11\t2\tENV\t256318\n",
  "Metarhizium anisopliae\t\t1\t4\tPLN\t5530\n",
  "Microbotryum violaceum\t\t1\t4\tPLN\t5272\n",
  "Microcebus murinus\tgray mouse lemur\t1\t2\tPRI\t30608\n",
  "Mimulus guttatus\tspotted monkey flower\t1\t1\tPLN\t4155\n",
  "Mimulus lewisii\t\t1\t1\tPLN\t69919\n",
  "Misgurnus anguillicaudatus\toriental weatherfish\t1\t2\tVRT\t75329\n",
  "Molgula tectiformis\t\t1\t13\tINV\t30286\n",
  "Monacrosporium haptotylum\t\t1\t4\tPLN\t74495\n",
  "Monodelphis domestica\tgray short-tailed opossum\t1\t2\tMAM\t13616\n",
  "Monosiga ovata\t\t1\t1\tINV\t81526\n",
  "mouse gut metagenome\t\t11\t2\tENV\t410661\n",
  "Mus musculus\thouse mouse\t1\t2\tROD\t10090\n",
  "Mus musculus domesticus\twestern European house mouse\t1\t2\tROD\t10092\n",
  "Mus musculus molossinus\tJapanese wild mouse\t1\t2\tROD\t57486\n",
  "Musa acuminata\tdessert banana\t1\t1\tPLN\t4641\n",
  "Mycobacterium tuberculosis\t\t11\t0\tBCT\t1773\n",
  "Mycosphaerella graminicola\t\t1\t4\tPLN\t54734\n",
  "Myotis lucifugus\tlittle brown bat\t1\t2\tMAM\t59463\n",
  "Mytilus californianus\tCalifornia mussel\t1\t5\tINV\t6549\n",
  "Mytilus galloprovincialis\tMediterranean mussel\t1\t5\tINV\t29158\n",
  "Myzus persicae\tgreen peach aphid\t1\t5\tINV\t13164\n",
  "Nakaseomyces delphensis\t\t1\t3\tPLN\t51657\n",
  "Nasonia giraulti\t\t1\t5\tINV\t7426\n",
  "Nasonia vitripennis\tjewel wasp\t1\t5\tINV\t7425\n",
  "Necator americanus\t\t1\t5\tINV\t51031\n",
  "Neisseria gonorrhoeae\t\t11\t0\tBCT\t485\n",
  "Neisseria meningitidis\t\t11\t0\tBCT\t487\n",
  "Nelumbo nucifera\t\t1\t1\tPLN\t4432\n",
  "Nematostella vectensis\tstarlet sea anemone\t1\t4\tINV\t45351\n",
  "Neosartorya fischeri NRRL 181\t\t1\t4\tPLN\t331117\n",
  "Neospora caninum\t\t1\t4\tINV\t29176\n",
  "Neovison vison\t\t1\t2\tMAM\t452646\n",
  "Neurospora crassa\t\t1\t4\tPLN\t5141\n",
  "Neurospora crassa OR74A\t\t1\t4\tPLN\t367110\n",
  "Newcastle disease virus\t\t1\t0\tVRL\t11176\n",
  "Nicotiana benthamiana\t\t1\t1\tPLN\t4100\n",
  "Nicotiana sylvestris\twood tobacco\t1\t1\tPLN\t4096\n",
  "Nicotiana tabacum\tcommon tobacco\t1\t1\tPLN\t4097\n",
  "Nippostrongylus brasiliensis\t\t1\t5\tINV\t27835\n",
  "Nuphar advena\t\t1\t1\tPLN\t77108\n",
  "Ochotona princeps\tAmerican pika\t1\t2\tMAM\t9978\n",
  "Ocimum basilicum\tsweet basil\t1\t1\tPLN\t39350\n",
  "Onchocerca volvulus\t\t1\t5\tINV\t6282\n",
  "Oncometopia nigricans\t\t1\t5\tINV\t266772\n",
  "Oncorhynchus mykiss\trainbow trout\t1\t2\tVRT\t8022\n",
  "Oncorhynchus nerka\tsockeye salmon\t1\t2\tVRT\t8023\n",
  "Oncorhynchus tshawytscha\tChinook salmon\t1\t2\tVRT\t74940\n",
  "Ophiostoma piliferum\t\t1\t4\tPLN\t38032\n",
  "Opisthorchis viverrini\t\t1\t9\tINV\t6198\n",
  "Oreochromis niloticus\tNile tilapia\t1\t2\tVRT\t8128\n",
  "Ornithorhynchus anatinus\tplatypus\t1\t2\tMAM\t9258\n",
  "Oryctolagus cuniculus\trabbit\t1\t2\tMAM\t9986\n",
  "Oryza alta\t\t1\t1\tPLN\t52545\n",
  "Oryza australiensis\t\t1\t1\tPLN\t4532\n",
  "Oryza brachyantha\t\t1\t1\tPLN\t4533\n",
  "Oryza coarctata\t\t1\t1\tPLN\t77588\n",
  "Oryza glaberrima\tAfrican rice\t1\t1\tPLN\t4538\n",
  "Oryza granulata\t\t1\t1\tPLN\t110450\n",
  "Oryza minuta\t\t1\t1\tPLN\t63629\n",
  "Oryza nivara\t\t1\t1\tPLN\t4536\n",
  "Oryza officinalis\t\t1\t1\tPLN\t4535\n",
  "Oryza punctata\t\t1\t1\tPLN\t4537\n",
  "Oryza ridleyi\t\t1\t1\tPLN\t83308\n",
  "Oryza rufipogon\t\t1\t1\tPLN\t4529\n",
  "Oryza sativa\trice\t1\t1\tPLN\t4530\n",
  "Oryza sativa Indica Group\t\t1\t1\tPLN\t39946\n",
  "Oryza sativa Japonica Group\t\t1\t1\tPLN\t39947\n",
  "Oryzias latipes\tJapanese medaka\t1\t2\tVRT\t8090\n",
  "Oscarella carmela\t\t1\t4\tINV\t386100\n",
  "Osmerus mordax\trainbow smelt\t1\t2\tVRT\t8014\n",
  "Ostertagia ostertagi\t\t1\t5\tINV\t6317\n",
  "Ostreococcus lucimarinus CCE9901\t\t1\t1\tPLN\t436017\n",
  "Otolemur garnettii\tsmall-eared galago\t1\t2\tPRI\t30611\n",
  "Ovis aries\tsheep\t1\t2\tMAM\t9940\n",
  "Pan paniscus\tpygmy chimpanzee\t1\t2\tPRI\t9597\n",
  "Pan troglodytes\tchimpanzee\t1\t2\tPRI\t9598\n",
  "Pan troglodytes troglodytes\t\t1\t2\tPRI\t37011\n",
  "Pan troglodytes verus\t\t1\t2\tPRI\t37012\n",
  "Panax ginseng\t\t1\t1\tPLN\t4054\n",
  "Panicum virgatum\tswitchgrass\t1\t1\tPLN\t38727\n",
  "Paracentrotus lividus\tcommon urchin\t1\t9\tINV\t7656\n",
  "Paracoccidioides brasiliensis\t\t1\t4\tPLN\t121759\n",
  "Paralabidochromis chilotes\t\t1\t2\tVRT\t77306\n",
  "Paralichthys olivaceus\tbastard halibut\t1\t2\tVRT\t8255\n",
  "Paramecium tetraurelia\t\t6\t4\tINV\t5888\n",
  "Parastrongyloides trichosuri\t\t1\t5\tINV\t131310\n",
  "Pasteuria penetrans\t\t11\t0\tBCT\t86005\n",
  "Patiria pectinifera\t\t1\t9\tINV\t7594\n",
  "Pavlova lutheri\t\t1\t1\tPLN\t2832\n",
  "Paxillus involutus\t\t1\t4\tPLN\t71150\n",
  "Pediculus humanus capitis\thuman head louse\t1\t5\tINV\t121226\n",
  "Pediculus humanus corporis\thuman body louse\t1\t5\tINV\t121224\n",
  "Penaeus monodon\tblack tiger shrimp\t1\t5\tINV\t6687\n",
  "Penicillium marneffei\t\t1\t4\tPLN\t37727\n",
  "Pennisetum glaucum\t\t1\t1\tPLN\t4543\n",
  "Perkinsus marinus ATCC 50983\t\t1\t4\tINV\t423536\n",
  "Peromyscus maniculatus bairdii\tprairie deer mouse\t1\t2\tROD\t230844\n",
  "Persea americana\tavocado\t1\t1\tPLN\t3435\n",
  "Petromyzon marinus\tsea lamprey\t1\t2\tVRT\t7757\n",
  "Phaeodactylum tricornutum\t\t1\t1\tPLN\t2850\n",
  "Phaeosphaeria nodorum SN15\t\t1\t4\tPLN\t321614\n",
  "Phakopsora pachyrhizi\t\t1\t4\tPLN\t170000\n",
  "Phalaenopsis equestris\t\t1\t1\tPLN\t78828\n",
  "Phalaenopsis violacea\t\t1\t1\tPLN\t86509\n",
  "Phanerochaete chrysosporium\t\t1\t4\tPLN\t5306\n",
  "Phaseolus coccineus\t\t1\t1\tPLN\t3886\n",
  "Phaseolus vulgaris\t\t1\t1\tPLN\t3885\n",
  "Phlebotomus papatasi\t\t1\t5\tINV\t29031\n",
  "Photorhabdus luminescens\t\t11\t0\tBCT\t29488\n",
  "Physarum polycephalum\tslime mold\t1\t1\tINV\t5791\n",
  "Physcomitrella patens\t\t1\t1\tPLN\t3218\n",
  "Physcomitrella patens subsp. patens\t\t1\t1\tPLN\t145481\n",
  "Phytophthora brassicae\t\t1\t1\tPLN\t187813\n",
  "Phytophthora infestans\tpotato late blight agent\t1\t1\tPLN\t4787\n",
  "Phytophthora infestans T30-4\t\t1\t1\tPLN\t403677\n",
  "Phytophthora parasitica\t\t1\t1\tPLN\t4792\n",
  "Phytophthora ramorum\tSudden oak death agent\t1\t1\tPLN\t164328\n",
  "Phytophthora sojae\t\t1\t1\tPLN\t67593\n",
  "Picea abies\tNorway spruce\t1\t1\tPLN\t3329\n",
  "Picea glauca\twhite spruce\t1\t1\tPLN\t3330\n",
  "Picea sitchensis\tSitka spruce\t1\t1\tPLN\t3332\n",
  "Pichia angusta\t\t1\t3\tPLN\t4905\n",
  "Pichia farinosa\t\t1\t3\tPLN\t4920\n",
  "Pichia guilliermondii ATCC 6260\t\t1\t3\tPLN\t294746\n",
  "Pichia stipitis CBS 6054\t\t12\t3\tPLN\t322104\n",
  "Pimephales promelas\t\t1\t2\tVRT\t90988\n",
  "Pinus pinaster\t\t1\t1\tPLN\t71647\n",
  "Pinus taeda\tloblolly pine\t1\t1\tPLN\t3352\n",
  "Pisum sativum\tpea\t1\t1\tPLN\t3888\n",
  "Plasmodium berghei\t\t1\t4\tINV\t5821\n",
  "Plasmodium berghei strain ANKA\t\t1\t4\tINV\t5823\n",
  "Plasmodium chabaudi\t\t1\t4\tINV\t5825\n",
  "Plasmodium chabaudi chabaudi\t\t1\t4\tINV\t31271\n",
  "Plasmodium falciparum\tmalaria parasite P. falciparum\t1\t4\tINV\t5833\n",
  "Plasmodium falciparum 3D7\t\t1\t4\tINV\t36329\n",
  "Plasmodium falciparum Dd2\t\t1\t4\tINV\t57267\n",
  "Plasmodium falciparum HB3\t\t1\t4\tINV\t137071\n",
  "Plasmodium vivax\tmalaria parasite P. vivax\t1\t4\tINV\t5855\n",
  "Plasmodium yoelii\t\t1\t4\tINV\t5861\n",
  "Plasmodium yoelii yoelii\t\t1\t4\tINV\t73239\n",
  "Plasmodium yoelii yoelii str. 17XNL\t\t1\t4\tINV\t352914\n",
  "Platichthys flesus\tEuropean flounder\t1\t2\tVRT\t8260\n",
  "Plodia interpunctella\tIndianmeal moth\t1\t5\tINV\t58824\n",
  "Pneumocystis carinii\t\t1\t4\tPLN\t4754\n",
  "Podocoryne carnea\t\t1\t4\tINV\t6096\n",
  "Poecilia reticulata\tguppy\t1\t2\tVRT\t8081\n",
  "Polygonatum sibiricum\t\t1\t1\tPLN\t261423\n",
  "Polysphondylium pallidum\t\t1\t1\tINV\t13642\n",
  "Polytomella parva\t\t1\t1\tPLN\t51329\n",
  "Poncirus trifoliata\t\t1\t1\tPLN\t37690\n",
  "Pongo pygmaeus\torangutan\t1\t2\tPRI\t9600\n",
  "Populus deltoides\t\t1\t1\tPLN\t3696\n",
  "Populus euphratica\t\t1\t1\tPLN\t75702\n",
  "Populus nigra\t\t1\t1\tPLN\t3691\n",
  "Populus tremula\t\t1\t1\tPLN\t113636\n",
  "Populus tremuloides\tquaking aspen\t1\t1\tPLN\t3693\n",
  "Populus trichocarpa\t\t1\t1\tPLN\t3694\n",
  "Porcine respiratory and reproductive syndrome virus\t\t1\t0\tVRL\t28344\n",
  "Porphyra haitanensis\t\t1\t4\tPLN\t76159\n",
  "Porphyra yezoensis\t\t1\t4\tPLN\t2788\n",
  "Pratylenchus vulnus\t\t1\t5\tINV\t45931\n",
  "Pristionchus pacificus\t\t1\t5\tINV\t54126\n",
  "Prototheca wickerhamii\t\t1\t1\tPLN\t3111\n",
  "Prunus armeniaca\tapricot\t1\t1\tPLN\t36596\n",
  "Prunus dulcis\talmond\t1\t1\tPLN\t3755\n",
  "Prunus persica\tpeach\t1\t1\tPLN\t3760\n",
  "Prymnesium parvum\t\t1\t1\tPLN\t97485\n",
  "Pseudomonas aeruginosa\t\t11\t0\tBCT\t287\n",
  "Pseudotsuga menziesii\tDouglas-fir\t1\t1\tPLN\t3357\n",
  "Pseudotsuga menziesii var. menziesii\t\t1\t1\tPLN\t278161\n",
  "Psychroflexus torquis ATCC 700755\t\t11\t0\tBCT\t313595\n",
  "Puccinellia tenuiflora\t\t1\t1\tPLN\t240906\n",
  "Puccinia graminis f. sp. tritici CRL 75-36-700-3\t\t1\t4\tPLN\t418459\n",
  "Puccinia triticina\t\t1\t4\tPLN\t208348\n",
  "Pythium ultimum DAOM BR144\t\t1\t1\tPLN\t431595\n",
  "Rabies virus\t\t1\t0\tVRL\t11292\n",
  "Raphanus raphanistrum subsp. landra\t\t1\t1\tPLN\t328428\n",
  "Raphanus raphanistrum subsp. maritimus\t\t1\t1\tPLN\t457179\n",
  "Raphanus raphanistrum subsp. raphanistrum\t\t1\t1\tPLN\t109997\n",
  "Raphanus sativus\tradish\t1\t1\tPLN\t3726\n",
  "Raphanus sativus var. oleiformis\t\t1\t1\tPLN\t463157\n",
  "Rattus norvegicus\tNorway rat\t1\t2\tROD\t10116\n",
  "Reclinomonas americana\t\t1\t1\tINV\t48483\n",
  "Rhipicephalus appendiculatus\t\t1\t5\tINV\t34631\n",
  "Rhipicephalus microplus\tsouthern cattle tick\t1\t5\tINV\t6941\n",
  "Rhizopus oryzae\t\t1\t1\tPLN\t64495\n",
  "Rhodnius prolixus\t\t1\t5\tINV\t13249\n",
  "Rhynchosciara americana\t\t1\t5\tINV\t7186\n",
  "Ricinus communis\tcastor bean\t1\t1\tPLN\t3988\n",
  "Robinia pseudoacacia\t\t1\t1\tPLN\t35938\n",
  "Rosa hybrid cultivar\t\t1\t1\tPLN\t128735\n",
  "Rutilus rutilus\troach minnow\t1\t2\tVRT\t48668\n",
  "Saccharomyces castellii\t\t1\t3\tPLN\t27288\n",
  "Saccharomyces cerevisiae\tbaker's yeast\t1\t3\tPLN\t4932\n",
  "Saccharomyces mikatae IFO 1815\t\t1\t3\tPLN\t226126\n",
  "Saccharomyces pastorianus\t\t1\t3\tPLN\t27292\n",
  "Saccharomyces servazzii\t\t1\t3\tPLN\t27293\n",
  "Saccharomyces uvarum\t\t1\t3\tPLN\t230603\n",
  "Saccharum hybrid cultivar\t\t1\t1\tPLN\t128810\n",
  "Saccharum officinarum\t\t1\t1\tPLN\t4547\n",
  "Saccoglossus kowalevskii\t\t1\t9\tINV\t10224\n",
  "Saitoella complicata\t\t1\t4\tPLN\t5606\n",
  "Salmo salar\tAtlantic salmon\t1\t2\tVRT\t8030\n",
  "Salvelinus fontinalis\tbrook trout\t1\t2\tVRT\t8038\n",
  "Salvia miltiorrhiza\t\t1\t1\tPLN\t226208\n",
  "Sarcocystis falcatula\t\t1\t4\tINV\t32593\n",
  "Sarcocystis neurona\t\t1\t4\tINV\t42890\n",
  "Saruma henryi\t\t1\t1\tPLN\t13258\n",
  "Sawyeria marylandensis\t\t1\t1\tINV\t194530\n",
  "Scenedesmus obliquus\t\t1\t22\tPLN\t3088\n",
  "Schistosoma japonicum\t\t1\t9\tINV\t6182\n",
  "Schistosoma mansoni\t\t1\t9\tINV\t6183\n",
  "Schizosaccharomyces pombe\tfission yeast\t1\t4\tPLN\t4896\n",
  "Schizosaccharomyces pombe 972h-\t\t1\t4\tPLN\t284812\n",
  "Schmidtea mediterranea\t\t1\t9\tINV\t79327\n",
  "Sclerotinia sclerotiorum\t\t1\t4\tPLN\t5180\n",
  "Sclerotinia sclerotiorum 1980\t\t1\t4\tPLN\t325569\n",
  "Sebastes rastrelliger\tgrass rockfish\t1\t2\tVRT\t72095\n",
  "Secale cereale\trye\t1\t1\tPLN\t4550\n",
  "Seculamonas ecuadoriensis\t\t1\t1\tINV\t221724\n",
  "Selaginella lepidophylla\t\t1\t1\tPLN\t59777\n",
  "Sesamum indicum\tsesame\t1\t1\tPLN\t4182\n",
  "Simian immunodeficiency virus\t\t1\t0\tVRL\t11723\n",
  "Solanum chacoense\tChaco potato\t1\t1\tPLN\t4108\n",
  "Solanum habrochaites\t\t1\t1\tPLN\t62890\n",
  "Solanum lycopersicum\t\t1\t1\tPLN\t4081\n",
  "Solanum pennellii\t\t1\t1\tPLN\t28526\n",
  "Solanum tuberosum\tpotato\t1\t1\tPLN\t4113\n",
  "Solenopsis invicta\tred fire ant\t1\t5\tINV\t13686\n",
  "Sorex araneus\tEuropean shrew\t1\t2\tMAM\t42254\n",
  "Sorghum bicolor\tsorghum\t1\t1\tPLN\t4558\n",
  "Sorghum propinquum\t\t1\t1\tPLN\t132711\n",
  "Sparus aurata\tgilthead seabream\t1\t2\tVRT\t8175\n",
  "Spermophilus lateralis\tgolden-mantled ground squirrel\t1\t2\tROD\t76772\n",
  "Spermophilus tridecemlineatus\tthirteen-lined ground squirrel\t1\t2\tROD\t43179\n",
  "Sphaeroforma arctica\t\t1\t1\tINV\t72019\n",
  "Spironucleus barkhanus\t\t6\t0\tINV\t103874\n",
  "Spizellomyces punctatus\t\t1\t16\tPLN\t109760\n",
  "Spodoptera frugiperda\tfall armyworm\t1\t5\tINV\t7108\n",
  "Squalus acanthias\tspiny dogfish\t1\t2\tVRT\t7797\n",
  "Stachyamoeba lipophora\t\t1\t1\tINV\t463046\n",
  "Staphylococcus aureus\t\t11\t0\tBCT\t1280\n",
  "Staphylococcus epidermidis\t\t11\t0\tBCT\t1282\n",
  "Sterkiella histriomuscorum\t\t6\t4\tINV\t94289\n",
  "Stevia rebaudiana\t\t1\t1\tPLN\t55670\n",
  "Streblomastix strix\t\t6\t1\tINV\t222440\n",
  "Streptococcus agalactiae\t\t11\t0\tBCT\t1311\n",
  "Streptococcus pneumoniae\t\t11\t0\tBCT\t1313\n",
  "Streptococcus pyogenes\t\t11\t0\tBCT\t1314\n",
  "Strongylocentrotus purpuratus\t\t1\t9\tINV\t7668\n",
  "Strongyloides ratti\t\t1\t5\tINV\t34506\n",
  "Strongyloides stercoralis\t\t1\t5\tINV\t6248\n",
  "Suidasia medanensis\t\t1\t5\tINV\t223625\n",
  "Sus scrofa\tpig\t1\t2\tMAM\t9823\n",
  "Taenia solium\tpork tapeworm\t1\t9\tINV\t6204\n",
  "Taeniopygia guttata\t\t1\t2\tVRT\t59729\n",
  "Taiwania cryptomerioides\t\t1\t1\tPLN\t50187\n",
  "Takifugu rubripes\ttorafugu\t1\t2\tVRT\t31033\n",
  "Tamarix androssowii\t\t1\t1\tPLN\t189785\n",
  "Tamarix hispida\t\t1\t1\tPLN\t189793\n",
  "Taphrina deformans\tpeach leaf curl fungus\t1\t4\tPLN\t5011\n",
  "Taraxacum kok-saghyz\t\t1\t1\tPLN\t333970\n",
  "Taraxacum officinale\t\t1\t1\tPLN\t50225\n",
  "Teladorsagia circumcincta\t\t1\t5\tINV\t45464\n",
  "Tetrahymena thermophila\t\t6\t4\tINV\t5911\n",
  "Tetrahymena thermophila SB210\t\t6\t4\tINV\t312017\n",
  "Tetraodon nigroviridis\t\t1\t2\tVRT\t99883\n",
  "Thalassiosira pseudonana CCMP1335\t\t1\t4\tPLN\t296543\n",
  "Theileria annulata\t\t1\t4\tINV\t5874\n",
  "Theileria annulata strain Ankara\t\t1\t4\tINV\t353154\n",
  "Theileria parva\t\t1\t4\tINV\t5875\n",
  "Theileria parva strain Muguga\t\t1\t4\tINV\t333668\n",
  "Thellungiella halophila\t\t1\t1\tPLN\t98038\n",
  "Thellungiella salsuginea\t\t1\t1\tPLN\t72664\n",
  "Theobroma cacao\tcacao\t1\t1\tPLN\t3641\n",
  "Thermomyces lanuginosus\t\t1\t4\tPLN\t5541\n",
  "Thlaspi caerulescens\t\t1\t1\tPLN\t107243\n",
  "Thunnus thynnus\tbluefin tuna\t1\t2\tVRT\t8237\n",
  "Thymallus arcticus\tArctic grayling\t1\t2\tVRT\t70285\n",
  "Torque teno virus\t\t1\t0\tVRL\t68887\n",
  "Tortula ruralis\t\t1\t1\tPLN\t38588\n",
  "Toxocara canis\t\t1\t5\tINV\t6265\n",
  "Toxoplasma gondii\t\t1\t4\tINV\t5811\n",
  "Toxoplasma gondii RH\t\t1\t4\tINV\t383379\n",
  "Toxoptera citricida\tbrown citrus aphid\t1\t5\tINV\t223852\n",
  "Trametes versicolor\t\t1\t4\tPLN\t5325\n",
  "Tribolium castaneum\tred flour beetle\t1\t5\tINV\t7070\n",
  "Trichinella spiralis\t\t1\t5\tINV\t6334\n",
  "Trichoderma reesei QM6a\t\t1\t4\tPLN\t431241\n",
  "Trichomonas vaginalis\t\t1\t0\tINV\t5722\n",
  "Trichomonas vaginalis G3\t\t1\t0\tINV\t412133\n",
  "Trichophyton rubrum\t\t1\t4\tPLN\t5551\n",
  "Trichosurus vulpecula\tsilver-gray brushtail possum\t1\t2\tMAM\t9337\n",
  "Trichuris muris\t\t1\t5\tINV\t70415\n",
  "Trichuris vulpis\t\t1\t5\tINV\t219738\n",
  "Trifolium pratense\t\t1\t1\tPLN\t57577\n",
  "Trifolium repens\twhite clover\t1\t1\tPLN\t3899\n",
  "Trimastix pyriformis\t\t1\t1\tINV\t342808\n",
  "Triphysaria versicolor\t\t1\t1\tPLN\t64093\n",
  "Triticum aestivum\tbread wheat\t1\t1\tPLN\t4565\n",
  "Triticum monococcum\t\t1\t1\tPLN\t4568\n",
  "Triticum turgidum subsp. durum\tdurum wheat\t1\t1\tPLN\t4567\n",
  "Tritrichomonas foetus\t\t1\t0\tINV\t5724\n",
  "Trypanosoma brucei\t\t1\t4\tINV\t5691\n",
  "Trypanosoma brucei rhodesiense\t\t1\t4\tINV\t31286\n",
  "Trypanosoma brucei TREU927\t\t1\t4\tINV\t185431\n",
  "Trypanosoma cruzi\t\t1\t4\tINV\t5693\n",
  "Trypanosoma cruzi strain CL Brener\t\t1\t4\tINV\t353153\n",
  "Tuber borchii\twhitish truffle\t1\t4\tPLN\t42251\n",
  "Tupaia belangeri\tnorthern tree shrew\t1\t2\tMAM\t37347\n",
  "Tursiops truncatus\tbottlenosed dolphin\t1\t2\tMAM\t9739\n",
  "uncultured archaeon\t\t11\t0\tENV\t115547\n",
  "uncultured bacterium\t\t11\t0\tENV\t77133\n",
  "uncultured eukaryote\t\t1\t1\tENV\t100272\n",
  "uncultured fungus\t\t1\t4\tENV\t175245\n",
  "uncultured marine virus\t\t1\t0\tENV\t186617\n",
  "uncultured organism\t\t11\t2\tENV\t155900\n",
  "uncultured prokaryote\t\t11\t0\tENV\t198431\n",
  "Uromyces appendiculatus\t\t1\t4\tPLN\t5264\n",
  "Ustilago maydis\t\t1\t4\tPLN\t5270\n",
  "Ustilago maydis 521\t\t1\t4\tPLN\t237631\n",
  "Vaccinium corymbosum\t\t1\t1\tPLN\t69266\n",
  "Vanderwaltozyma polyspora DSM 70294\t\t1\t3\tPLN\t436907\n",
  "Verticillium dahliae\t\t1\t4\tPLN\t27337\n",
  "Vibrio parahaemolyticus AQ3810\t\t11\t0\tBCT\t419109\n",
  "Vigna unguiculata\tcowpea\t1\t1\tPLN\t3917\n",
  "Vitis hybrid cultivar\t\t1\t1\tPLN\t241073\n",
  "Vitis shuttleworthii\t\t1\t1\tPLN\t246827\n",
  "Vitis vinifera\t\t1\t1\tPLN\t29760\n",
  "Welwitschia mirabilis\t\t1\t1\tPLN\t3377\n",
  "Wuchereria bancrofti\t\t1\t5\tINV\t6293\n",
  "Xenopus laevis\tAfrican clawed frog\t1\t2\tVRT\t8355\n",
  "Xenopus tropicalis\twestern clawed frog\t1\t2\tVRT\t8364\n",
  "Xiphinema index\t\t1\t5\tINV\t46003\n",
  "Yarrowia lipolytica\t\t1\t3\tPLN\t4952\n",
  "Yarrowia lipolytica CLIB122\t\t1\t3\tPLN\t284591\n",
  "Zamia fischeri\t\t1\t1\tPLN\t34342\n",
  "Zantedeschia aethiopica\t\t1\t1\tPLN\t69721\n",
  "Zea mays\t\t1\t1\tPLN\t4577\n",
  "Zea mays subsp. mays\tmaize\t1\t1\tPLN\t381124\n",
  "Zea mays subsp. parviglumis\t\t1\t1\tPLN\t76912\n",
  "Zingiber officinale\t\t1\t1\tPLN\t94328\n",
  "Zinnia elegans\t\t1\t1\tPLN\t34245\n",
  "Zostera marina\t\t1\t1\tPLN\t29655\n",
  "Zygosaccharomyces rouxii\t\t1\t3\tPLN\t4956\n",
  NULL
};
#endif

extern Boolean LoadOrganismTable (void)

{
  Char     ch;
  FILE     *f;
  Boolean  failed;
  Char     first [16];
  Int2     idx;
  CharPtr  lst;
  Int8     len;
  Int2     num;
  CharPtr  ptr = NULL;
  Char     str [PATH_MAX];
  CharPtr  tmp;
  Int2     version;
#ifdef WIN_MAC
  CharPtr  p;
#endif
#if (defined(OS_DOS) || defined (OS_NT))
  CharPtr  p;
  CharPtr  q;
#endif

  orgTxtPtr = NULL;
  orgStrIdx = NULL;
  orgNum = 0;
  failed = TRUE;
  ProgramPath (str, sizeof (str));
  tmp = StringRChr (str, DIRDELIMCHR);
  if (tmp != NULL) {
    *tmp = '\0';
    FileBuildPath (str, NULL, "taxlist.txt");
    len = FileLength (str);
    f = FileOpen (str, "r");
    if (f == NULL) {
      if (GetAppParam ("NCBI", "NCBI", "DATA", "", str, sizeof (str))) {
        FileBuildPath (str, NULL, "taxlist.txt");
        len = FileLength (str);
        f = FileOpen (str, "r");
      }
    }
    if (f != NULL) {
      ptr = MemNew ((size_t) (len + 5));
      if (ptr != NULL) {
        FileRead (ptr, (size_t) len, 1, f);
#if (defined(OS_DOS) || defined (OS_NT))
        p = ptr;
        q = ptr;
        while (*p) {
          if (*p == '\r') {
            p++;
          } else {
            *q = *p;
            p++;
            q++;
          }
        }
        *q = '\0';
#endif
#ifdef WIN_MAC
        p = ptr;
        while (*p) {
          if (*p == '\r') {
            *p = '\n';
          }
          p++;
        }
#endif
      }
      FileClose (f);
    }
  }
  if (ptr == NULL) {
#ifndef WIN16
    idx = 0;
    len = 0;
    tmp = taxlistMemStrs [idx];
    while (tmp != NULL) {
      len += StringLen (tmp);
      idx++;
      tmp = taxlistMemStrs [idx];
    }
    ptr = MemNew ((size_t) (len + 5));
    if (ptr != NULL) {
      lst = ptr;
      idx = 0;
      tmp = taxlistMemStrs [idx];
      while (tmp != NULL) {
        lst = StringMove (lst, tmp);
        idx++;
        tmp = taxlistMemStrs [idx];
      }
    }
#endif
  }
  if (ptr != NULL) {
    orgTxtPtr = ptr;
    tmp = ptr;
    ch = *tmp;
    while (ch != '\0' && ch != '\n') {
      tmp++;
      ch = *tmp;
    }
    *tmp = '\0';
    StringNCpy_0 (first, ptr, sizeof (first) - 1);
    *tmp = ch;
    if (StrToInt (first, &version)) {
      tmp++;
      ptr = tmp;
    }
    num = 0;
    tmp = ptr;
    ch = *tmp;
    while (ch != '\0') {
      if (ch == '\n') {
        num++;
      }
      tmp++;
      ch = *tmp;
    }
    orgStrIdx = MemNew (sizeof (CharPtr) * (size_t) (num + 3));
    if (orgStrIdx != NULL) {
      idx = 0;
      tmp = ptr;
      ch = *tmp;
      orgStrIdx [idx] = tmp;
      while (ch != '\0') {
        if (ch == '\n') {
          idx++;
          tmp++;
          ch = *tmp;
          orgStrIdx [idx] = tmp;
        } else {
          tmp++;
          ch = *tmp;
        }
      }
      orgNum = num;
      failed = FALSE;
    }
  }
  if (failed) {
    orgTxtPtr = MemFree (orgTxtPtr);
    orgStrIdx = MemFree (orgStrIdx);
    return FALSE;
  } else {
    return TRUE;
  }
}

extern void FreeOrganismTable (void)

{
  orgTxtPtr = MemFree (orgTxtPtr);
  orgStrIdx = MemFree (orgStrIdx);
}

extern void SetupGeneticCodes (void)

{
  Char            ch;
  GeneticCodePtr  codes;
  GeneticCodePtr  gcp;
  Int2            i;
  Int4            id;
  Int2            j;
  Int2            index;
  Char            name [64];
  CharPtr         ptr;
  Char            str [256];
  ValNodePtr      tmp;

  numGeneticCodes = 0;
  for (i = 0; i < NUM_GENETIC_CODES; i++) {
    gcIndexToId [i] = 0;
    gcIdToIndex [i] = 1;
    gcNames [i] = NULL;
  }
  index = 1;
  codes = GeneticCodeTableLoad ();
  if (codes != NULL) {
    for (gcp = codes; gcp != NULL; gcp = gcp->next) {
      id = 0;
      str [0] = '\0';
      for (tmp = (ValNodePtr) gcp->data.ptrvalue; tmp != NULL; tmp = tmp->next) {
        switch (tmp->choice) {
          case 1 :
            if (StringLen (str) < 1) {
              StringNCpy_0 (str, (CharPtr) tmp->data.ptrvalue, sizeof (str));
              ptr = str;
              ch = *ptr;
              while (ch != '\0') {
                if (ch == '/') {
                  *ptr = '-';
                }
                ptr++;
                ch = *ptr;
              }
            }
            break;
          case 2 :
            id = tmp->data.intvalue;
            break;
          default :
            break;
        }
      }
      if (id != 7 && id != 8) {
        if (id > 0 /* && id < 30 */ ) {
          i = 0;
          if (StringLen (str + i) > 0 && index < NUM_GENETIC_CODES - 1) {
            ch = str [i];
            while (ch == ' ' || ch == ';') {
              i++;
              ch = str [i];
            }
            j = 0;
            ch = str [i + j];
            while (ch != '\0' && ch != ';') {
              name [j] = ch;
              j++;
              ch = str [i + j];
            }
            name [j] = '\0';
            i += j;
            index++;
            if (ch == ';') {
              StringCat (name, ", etc.");
            }
            gcIndexToId [index] = id;
            gcIdToIndex [id] = index;
            gcNames [index] = StringSave (name);
          }
        }
      }
    }
  }
  numGeneticCodes = index;
}

extern void FreeGeneticCodes (void)

{
  Int2  i;

  for (i = 0; i < NUM_GENETIC_CODES; i++) {
    gcNames [i] = MemFree (gcNames [i]);
  }
}

extern ValNodePtr GetGeneticCodeValNodeList (void)
{
  ValNodePtr gencodelist = NULL;
  Int4       index;
  
  for (index = 0; index <= numGeneticCodes; index++)
  {
    if (StringHasNoText (gcNames[index]))
    {
      continue;
    }
    ValNodeAddPointer (&gencodelist, gcIndexToId [index], StringSave (gcNames[index]));
  }
  return gencodelist;
}

static void CopyField (CharPtr str, size_t max, CharPtr source, Int2 col)

{
  Char     ch;
  size_t   count;
  CharPtr  ptr;

  if (str != NULL && max > 0 && source != NULL) {
    MemSet (str, 0, max);
      ptr = source;
      ch = *ptr;
      while (col > 1 && ch != '\n' && ch != '\0') {
        while (ch != '\t' && ch != '\n' && ch != '\0') {
          ptr++;
          ch = *ptr;
        }
        if (ch == '\t') {
          ptr++;
          ch = *ptr;
        }
        col--;
      }
      count = 0;
      ch = ptr [count];
      while (ch != '\t' && ch != '\n' && ch != '\0') {
        count++;
        ch = ptr [count];
      }
      max = MIN (max, count);
      StringNCpy (str, ptr, max); /* remains StringNCpy, not _0 */
  }
}

static void CopyStrFromTaxPtr (CharPtr str, size_t max, Int2 row, Int2 col)

{
  CharPtr  source;

  if (str != NULL && max > 0) {
    MemSet (str, 0, max);
    if (orgTxtPtr != NULL && orgStrIdx != NULL && row > 0) {
      source = orgStrIdx [row - 1];
      CopyField (str, max, source, col);
    }
  }
}

static Int2 FindTaxText (CharPtr text, Int2 num)

{
  Int2  compare;
  Int2  left;
  Int2  mid;
  Int2  right;
  Char  str [256];

  mid = 0;
  if (text != NULL && num > 0) {
    left = 1;
    right = num;
    while (left <= right) {
      mid = (left + right) / 2;
      CopyStrFromTaxPtr (str, sizeof (str) - 2, mid, 1);
      compare = StringICmp (text, str);
      if (compare <= 0) {
        right = mid - 1;
      }
      if (compare >= 0) {
        left = mid + 1;
      }
    }
    if (left <= right + 1) {
      CopyStrFromTaxPtr (str, sizeof (str) - 2, mid, 1);
      str [StringLen (text)] = '\0';
      compare = StringICmp (text, str);
      if (compare > 0) {
        mid++;
        if (mid <= num) {
          CopyStrFromTaxPtr (str, sizeof (str) - 2, mid, 1);
          str [StringLen (text)] = '\0';
          compare = StringICmp (text, str);
          if (compare > 0) {
            mid = 0;
          }
        }
      }
    }
  }
  return mid;
}

extern Boolean BioSourceDialogToGenBankDivision (DialoG d, CharPtr div, size_t maxsize)

{
  GenBioPagePtr  gbp;

  if (div == NULL || maxsize < 1) return FALSE;
  div [0] = '\0';
  gbp = (GenBioPagePtr) GetObjectExtra (d);
  if (gbp != NULL) {
    GetTitle (gbp->gbDiv, div, maxsize);
    if (! StringHasNoText (div)) return TRUE;
  }
  return FALSE;
}

static void PopulateGeneticCodePopup (PopuP gc)

{
  Int2  i;

   if (gc != NULL) {
    PopupItem (gc, " ");
    for (i = 1; i <= numGeneticCodes; i++) {
      PopupItem (gc, gcNames [i]);
    }
  }
}

static void ChangeGencodePopups (GenBioPagePtr gbp)

{
  UIEnum  genome;
  ValNodePtr vnp;

  if (gbp != NULL) {
    if (gbp->simplecode != NULL) {
      vnp = DialogToPointer (gbp->genome);
      if (vnp != NULL) {
        genome = GenomeFromSrcLoc (vnp->choice);
        vnp = ValNodeFreeData (vnp);
        if (genome == 4 || genome == 5) {
          SafeSetValue (gbp->simplecode, gcIdToIndex [gbp->mitoGC]);
        } else if (genome == GENOME_chloroplast ||
                   genome == GENOME_chromoplast ||
                   genome == GENOME_plastid ||
                   genome == GENOME_cyanelle ||
                   genome == GENOME_apicoplast ||
                   genome == GENOME_leucoplast ||
                   genome == GENOME_proplastid) {
          SafeSetValue (gbp->simplecode, gcIdToIndex [11]);
        } else {
          SafeSetValue (gbp->simplecode, gcIdToIndex [gbp->nuclGC]);
        }
      }
      SafeSetValue (gbp->gcode, gcIdToIndex [gbp->nuclGC]);
      SafeSetValue (gbp->mgcode, gcIdToIndex [gbp->mitoGC]);
    } else {
      SafeSetValue (gbp->gcode, gcIdToIndex [gbp->nuclGC]);
      SafeSetValue (gbp->mgcode, gcIdToIndex [gbp->mitoGC]);
    }
  }
}

static void SetCodes (GenBioPagePtr gbp, Int2 row, Boolean changepopups)

{
  Char  str [256];

  if (gbp != NULL && row > 0) {
    CopyStrFromTaxPtr (str, sizeof (str) - 2, row, 3);
    StrToInt (str, &gbp->nuclGC);
    CopyStrFromTaxPtr (str, sizeof (str) - 2, row, 4);
    StrToInt (str, &gbp->mitoGC);
    CopyStrFromTaxPtr (str, sizeof (str) - 2, row, 5);
    SafeSetTitle (gbp->gbDiv, str);
    CopyStrFromTaxPtr (str, sizeof (str) - 2, row, 6);
    StrToLong (str, &gbp->taxID);
    if (changepopups) {
      ChangeGencodePopups (gbp);
    }
  }
}

extern void SetGenome (PopuP p);
void SetGenome (PopuP p)

{
  GenBioPagePtr  gbp;

  gbp = (GenBioPagePtr) GetObjectExtra (p);
  if (gbp != NULL) {
    ChangeGencodePopups (gbp);
  }
}

static void AutoScrollTax (GenBioPagePtr gbp, TexT t, Boolean isSciName,
                           Boolean setText, Boolean changepopups)

{
  Int2  num;
  Int2  oldOrg;
  Int2  row;
  Char  str [256];
  Char  txt [256];

  if (gbp != NULL && t != NULL) {
    gbp->nuclGC = 0;
    gbp->mitoGC = 0;
    gbp->taxID = 0;
    SafeSetTitle (gbp->gbDiv, "");
    GetTitle (t, str, sizeof (str) - 2);
    if (str [0] == '\0') {
      if (isSciName) {
        gbp->typedSciName = FALSE;
      } else {
        gbp->typedComName = FALSE;
      }
    }
    num = 0;
    oldOrg = gbp->selectedOrg;
    GetItemParams (gbp->orglist, 1, NULL, &num, NULL, NULL, NULL);
    if (num > 0) {
      row = FindTaxText (str, num);
      if (row > 0 && row <= num) {
        if (! RowIsVisible (gbp->orglist, 1, row, NULL, NULL)) {
          SetOffset (gbp->orglist, 0, row - 1);
        }
        CopyStrFromTaxPtr (txt, sizeof (txt) - 2, row, 1);
        if (StringICmp (txt, str) != 0) {
          if (setText) {
            if (isSciName) {
              if (! gbp->typedComName) {
                SafeSetTitle (gbp->commonName, "");
              }
            } else {
              if (! gbp->typedSciName) {
                /*
                SafeSetTitle (gbp->taxName, "");
                */
              }
            }
          }
          gbp->selectedOrg = 0;
          InvalDocRows (gbp->orglist, 1, oldOrg, oldOrg);
          if (changepopups) {
            ChangeGencodePopups (gbp);
          }
          return;
        }
        if (isSciName) {
          CopyStrFromTaxPtr (str, sizeof (str) - 2, row, 2);
          if (! gbp->typedComName) {
            if (setText) {
              SafeSetTitle (gbp->commonName, str);
            }
          } else {
            GetTitle (gbp->commonName, txt, sizeof (txt));
            if (StringICmp (txt, str) != 0) {
              if (changepopups) {
                ChangeGencodePopups (gbp);
              }
              return;
            }
          }
        } else {
          if (changepopups) {
            ChangeGencodePopups (gbp);
          }
          return;
        }
        SetCodes (gbp, row, changepopups);
        if (oldOrg != row) {
          gbp->selectedOrg = 0;
          InvalDocRows (gbp->orglist, 1, oldOrg, oldOrg);
          gbp->selectedOrg = row;
          InvalDocRows (gbp->orglist, 1, row, row);
        }
      }
    }
  }
}

static void FreeGenBioOrgNameData (GenBioPagePtr gbp);

static void TaxNameText (TexT t)

{
  DbtagPtr         dbt;
  GenBioPagePtr    gbp;
  ValNodePtr       head;
  ValNodePtr       nextvnp;
  ValNodePtr PNTR  prevvnp;
  ValNodePtr       vnp;

  gbp = (GenBioPagePtr) GetObjectExtra (t);
  if (gbp != NULL) {
    gbp->typedSciName = TRUE;
    AutoScrollTax (gbp, t, TRUE, TRUE, TRUE);
    SafeSetTitle (gbp->lineage, "");
    SafeSetTitle (gbp->gbDiv, "");
    head = DialogToPointer (gbp->db);
    if (head != NULL) {
      prevvnp = &head;
      vnp = head;
      while (vnp != NULL) {
        nextvnp = vnp->next;
        dbt = (DbtagPtr) vnp->data.ptrvalue;
        if (dbt != NULL) {
          if (StringICmp (dbt->db, "taxon") == 0) {
            *prevvnp = vnp->next;
            vnp->next = NULL;
            DbtagFree (dbt);
            ValNodeFree (vnp);
            FreeGenBioOrgNameData (gbp);
          } else {
            prevvnp = (ValNodePtr PNTR) &(vnp->next);
          }
        }
        vnp = nextvnp;
      }
      PointerToDialog (gbp->db, head);
    }
  }
}

static void CommonNameText (TexT t)

{
  GenBioPagePtr  gbp;

  gbp = (GenBioPagePtr) GetObjectExtra (t);
  if (gbp != NULL) {
    /*
    gbp->typedComName = TRUE;
    AutoScrollTax (gbp, t, FALSE, TRUE, TRUE);
    */
  }
}

static void ClickTaxName (DoC d, PoinT pt)

{
  GenBioPagePtr  gbp;
  Int2           item;
  Int2           oldOrg;
  Int2           row;

  gbp = (GenBioPagePtr) GetObjectExtra (d);
  if (gbp != NULL) {
    MapDocPoint (d, pt, &item, &row, NULL, NULL);
    if (item > 0 && row > 0) {
      gbp->clickedOrg = row;
      oldOrg = gbp->selectedOrg;
      if (gbp->selectedOrg != row) {
        gbp->selectedOrg = 0;
        InvalDocRows (d, 1, oldOrg, oldOrg);
        gbp->selectedOrg = row;
        InvalDocRows (d, 1, row, row);
      }
    }
  }
}

static void ReleaseTaxName (DoC d, PoinT pt)

{
  GenBioPagePtr  gbp;
  Int2           item;
  Int2           row;
  Char           str [256];

  gbp = (GenBioPagePtr) GetObjectExtra (d);
  if (gbp != NULL) {
    gbp->nuclGC = 0;
    gbp->mitoGC = 0;
    gbp->taxID = 0;
    SafeSetTitle (gbp->gbDiv, "");
    MapDocPoint (d, pt, &item, &row, NULL, NULL);
    if (item > 0 && row > 0 && row == gbp->clickedOrg) {
      ResetClip ();
      if (row == gbp->clickedOrg) {
        CopyStrFromTaxPtr (str, sizeof (str) - 2, row, 1);
        SafeSetTitle (gbp->taxName, str);
        CopyStrFromTaxPtr (str, sizeof (str) - 2, row, 2);
        SafeSetTitle (gbp->commonName, str);
        Select (gbp->taxName);
        TaxNameText (gbp->taxName);
        SetCodes (gbp, row, TRUE);
      }
    }
  }
}

extern Boolean SetBioSourceDialogTaxName (DialoG d, CharPtr taxname)

{
  GenBioPagePtr gbp;
  Int2          num;
  Int2          oldOrg;
  Int2          row;
  Char          str [256];

  gbp = (GenBioPagePtr) GetObjectExtra (d);
  if (gbp == NULL) return FALSE;
  num = 0;
  oldOrg = gbp->selectedOrg;
  GetItemParams (gbp->orglist, 1, NULL, &num, NULL, NULL, NULL);
  if (num > 0) {
    row = FindTaxText (taxname, num);
    if (row > 0 && row <= num) {
      if (! RowIsVisible (gbp->orglist, 1, row, NULL, NULL)) {
        SetOffset (gbp->orglist, 0, row - 1);
      }
      CopyStrFromTaxPtr (str, sizeof (str) - 2, row, 1);
      if (StringICmp (str, taxname) != 0) {
        SafeSetTitle (gbp->commonName, "");
        SafeSetTitle (gbp->taxName, taxname);
        gbp->selectedOrg = 0;
        InvalDocRows (gbp->orglist, 1, oldOrg, oldOrg);
        ChangeGencodePopups (gbp);
        return TRUE;
      }
      CopyStrFromTaxPtr (str, sizeof (str) - 2, row, 1);
      SafeSetTitle (gbp->taxName, str);
      CopyStrFromTaxPtr (str, sizeof (str) - 2, row, 2);
      SafeSetTitle (gbp->commonName, str);
      Select (gbp->taxName);
      SetCodes (gbp, row, TRUE);
      gbp->selectedOrg = row;
      InvalDocRows (gbp->orglist, 1, oldOrg, oldOrg);
      InvalDocRows (gbp->orglist, 1, row, row);
      return TRUE;
    }
    else
    {
      SafeSetTitle (gbp->taxName, taxname);
      SafeSetTitle (gbp->commonName, "");
      gbp->selectedOrg = 0;
      InvalDocRows (gbp->orglist, 1, oldOrg, oldOrg);
      ChangeGencodePopups (gbp);
    }
  }
  return FALSE;
}

static Boolean HighlightTaxName (DoC d, Int2 item, Int2 row, Int2 col)

{
  GenBioPagePtr  gbp;

  gbp = (GenBioPagePtr) GetObjectExtra (d);
  if (gbp != NULL) {
    return (Boolean) (row == gbp->selectedOrg);
  } else {
    return FALSE;
  }
}

static ParData orgListPar = {FALSE, FALSE, FALSE, FALSE, FALSE, 0, 0};
static ColData orgListCol [] = {
  {0, 0, 80, 0, NULL, 'l', FALSE, FALSE, FALSE, FALSE, FALSE},
  {0, 0,  0, 0, NULL, 'l', FALSE, FALSE, FALSE, FALSE, FALSE},
  {0, 0,  0, 0, NULL, 'l', FALSE, FALSE, FALSE, FALSE, FALSE},
  {0, 0,  0, 0, NULL, 'l', FALSE, FALSE, FALSE, FALSE, FALSE},
  {0, 0,  0, 0, NULL, 'l', FALSE, FALSE, FALSE, FALSE, FALSE},
  {0, 0,  0, 0, NULL, 'l', FALSE, FALSE, FALSE, FALSE, FALSE},
  {0, 0,  0, 0, NULL, 'l', FALSE, FALSE, FALSE, FALSE, FALSE},
  {0, 0,  0, 0, NULL, 'l', FALSE, FALSE, FALSE, FALSE, TRUE}};

static CharPtr AllButFirstLinePrtProc (DoC d, Int2 item, Pointer ptr)

{
  Char     ch;
  CharPtr  tmp;

  if (ptr != NULL) {
    tmp = (CharPtr) ptr;
    ch = *tmp;
    while (ch != '\0' && ch != '\n') {
      tmp++;
      ch = *tmp;
    }
    if (ch != '\0') {
      tmp++;
      return StringSave (tmp);
    } else {
      return NULL;
    }
  } else {
    return NULL;
  }
}


static Char useGenomicText [] = "\
(Use 'Genomic' for a sequence encoded by a nuclear gene.)\n";

static Pointer MakeOrgNameDataCopy (OrgNamePtr onp, Uint1Ptr orgname_choiceP)
{
  OrgNamePtr onp_copy;
  Pointer    retval;

  if (orgname_choiceP != NULL) {
    *orgname_choiceP = 0;
  }

  if (onp == NULL) return NULL;

  onp_copy = (OrgNamePtr) AsnIoMemCopy (onp, (AsnReadFunc) OrgNameAsnRead, (AsnWriteFunc) OrgNameAsnWrite);
  if (onp_copy == NULL) return NULL;

  retval = onp_copy->data;
  if (orgname_choiceP != NULL) {
    *orgname_choiceP = onp_copy->choice;
  }
  onp_copy->data = NULL;
  onp_copy->choice = 0;
  OrgNameFree (onp_copy);
  return retval;
}

static void FreeGenBioOrgNameData (GenBioPagePtr gbp)
{
  OrgNamePtr     onp;

  if (gbp != NULL) {
    if (gbp->orgname_choice != 0) {
      onp = OrgNameNew ();
      if (onp != NULL) {
        onp->choice = gbp->orgname_choice;
        onp->data = gbp->orgname_data;
        OrgNameFree (onp);
        gbp->orgname_choice = 0;
        gbp->orgname_data = NULL;
      }
    }
  }
}


static void BioSourcePtrToGenBioPage (DialoG d, Pointer data)

{
  BioSourcePtr   biop;
  GenBioPagePtr  gbp;
  UIEnum         genome;
  OrgModPtr      mod;
  OrgNamePtr     onp;
  OrgRefPtr      orp;
  WindoW         tempPort;
  ValNode        vn;
  ValNodePtr     vnp;

  gbp = (GenBioPagePtr) GetObjectExtra (d);
  biop = (BioSourcePtr) data;
  if (gbp != NULL) {
    orp = NULL;
    onp = NULL;
    tempPort = SavePort (gbp->taxName);
    SafeSetTitle (gbp->taxName, "");
    SafeSetTitle (gbp->commonName, "");
    gbp->typedSciName = FALSE;
    gbp->typedComName = FALSE;
    gbp->selectedOrg = 0;
    gbp->clickedOrg = 0;
    gbp->nuclGC = 0;
    gbp->mitoGC = 0;
    gbp->taxID = 0;
    SafeSetTitle (gbp->gbDiv, "");
    SafeSetValue (gbp->genome, 2);
    SafeSetValue (gbp->origin, 1);
    SafeSetStatus (gbp->is_focus, FALSE);
    SafeSetValue (gbp->gcode, 1);
    SafeSetValue (gbp->mgcode, 1);
    SafeSetValue (gbp->simplecode, 1);
    SafeSetTitle (gbp->gbacr, "");
    SafeSetTitle (gbp->gbana, "");
    SafeSetTitle (gbp->gbsyn, "");
    SafeSetTitle (gbp->lineage, "");
    PointerToDialog (gbp->db, NULL);
    PointerToDialog (gbp->syn, NULL);
    PointerToDialog (gbp->mod, NULL);
    PointerToDialog (gbp->subsrc_val_dlg, NULL);
    PointerToDialog (gbp->orgmod_val_dlg, NULL);
    
    if (biop != NULL) {
      vn.choice = SrcLocFromGenome (biop->genome);
      vn.data.ptrvalue = NULL;
      vn.next = NULL;
      PointerToDialog (gbp->genome, &vn);
      SetEnumPopup (gbp->origin, biosource_origin_alist, (UIEnum) biop->origin);
      SafeSetStatus (gbp->is_focus, biop->is_focus);
      if (biop->is_focus) {
        SafeEnable (gbp->is_focus);
      }
      orp = biop->org;
      if (orp != NULL) {
        gbp->origTaxName = StringSave (orp->taxname);
        SafeSetTitle (gbp->taxName, orp->taxname);
        SafeSetTitle (gbp->commonName, orp->common);
        PointerToDialog (gbp->db, orp->db);
        PointerToDialog (gbp->syn, orp->syn);
        PointerToDialog (gbp->mod, orp->mod);
        onp = orp->orgname;
        if (onp != NULL) {
          /* store orgname data for unaltered retrieval later */
          gbp->orgname_choice = 0;
          gbp->orgname_data = MakeOrgNameDataCopy(onp, &gbp->orgname_choice);

          SafeSetTitle (gbp->lineage, onp->lineage);
          PointerToDialog (gbp->orgmod_val_dlg, onp->mod);
          mod = onp->mod;
          while (mod != NULL) {
            switch (mod->subtype) {
              case 32 :
                SetTitle (gbp->gbacr, mod->subname);
                break;
              case 33 :
                SetTitle (gbp->gbana, mod->subname);
                break;
              case 34 :
                SetTitle (gbp->gbsyn, mod->subname);
                break;
                break;
              default :
                break;
            }
            mod = mod->next;
          }
        }
      }

      PointerToDialog (gbp->subsrc_val_dlg, biop->subtype);
    }
    if (orp != NULL) {
      if (! TextHasNoText (gbp->taxName)) {
        AutoScrollTax (gbp, gbp->taxName, TRUE, FALSE, FALSE);
      /*
      } else if (! TextHasNoText (gbp->commonName)) {
        AutoScrollTax (gbp, gbp->commonName, FALSE, FALSE, FALSE);
      */
        SafeSetTitle (gbp->gbDiv, "");
        SafeSetValue (gbp->gcode, 1);
        SafeSetValue (gbp->mgcode, 1);
        SafeSetValue (gbp->simplecode, 1);
      }
    }
    if (onp != NULL) {
      SafeSetTitle (gbp->gbDiv, onp->div);
      if (onp->gcode > 0) {
        gbp->nuclGC = onp->gcode;
      }
      if (onp->mgcode > 0) {
        gbp->mitoGC = onp->mgcode;
      }
      if (gbp->simplecode != NULL) {
        vnp = DialogToPointer (gbp->genome);
        if (vnp != NULL) {
          genome = GenomeFromSrcLoc (vnp->choice);
          vnp = ValNodeFreeData (vnp);
          if (genome == 4 || genome == 5) {
            SafeSetValue (gbp->simplecode, gcIdToIndex [onp->mgcode]);
          } else if (genome == GENOME_chloroplast ||
                     genome == GENOME_chromoplast ||
                     genome == GENOME_plastid ||
                     genome == GENOME_cyanelle ||
                     genome == GENOME_apicoplast ||
                     genome == GENOME_leucoplast ||
                     genome == GENOME_proplastid) {
            SafeSetValue (gbp->simplecode, gcIdToIndex [11]);
          } else {
            SafeSetValue (gbp->simplecode, gcIdToIndex [onp->gcode]);
          }
        }
        SafeSetValue (gbp->gcode, gcIdToIndex [onp->gcode]);
        SafeSetValue (gbp->mgcode, gcIdToIndex [onp->mgcode]);
      } else {
        SafeSetValue (gbp->gcode, gcIdToIndex [onp->gcode]);
        SafeSetValue (gbp->mgcode, gcIdToIndex [onp->mgcode]);
      }
    }
    RestorePort (tempPort);
  }
}


static Pointer GenBioPageToBioSourcePtr (DialoG d)

{
  BioSourcePtr   biop;
  Char           buf [256];
  Char           ch;
  CharPtr        chptr;
  FILE           *f;
  GenBioPagePtr  gbp;
  UIEnum         genome;
  Boolean        goOn;
  OrgModPtr      mod;
  OrgModPtr      nextmod;
  OrgNamePtr     onp;
  OrgRefPtr      orp;
  OrgModPtr      PNTR prevmod;
  CharPtr        ptr;
  Char           str [PATH_MAX];
  Int4           taxID;
  OrgModPtr      tmpmod;
  SubSourcePtr   tmpssp;
  UIEnum         val;
  
  Int2           num; /* contains number of items in gbp->orglist */
  Int4           row; /* contains closest row to match in gbp->orglist */
  Char           txt [256]; /* holds tax name copied from gbp->orglist */
  ValNodePtr     vnp;

  biop = NULL;
  gbp = (GenBioPagePtr) GetObjectExtra (d);
  if (gbp != NULL) {
    biop = BioSourceNew ();
    if (biop != NULL) {

      vnp = DialogToPointer (gbp->genome);
      if (vnp != NULL) {
        biop->genome = (Uint1) GenomeFromSrcLoc (vnp->choice);
        vnp = ValNodeFreeData (vnp);
      }
      if (GetEnumPopup (gbp->origin, biosource_origin_alist, &val)) {
        biop->origin = (Uint1) val;
      }
      biop->is_focus = GetStatus (gbp->is_focus);
      orp = OrgRefNew ();
      biop->org = orp;
      if (orp != NULL) {
        orp->taxname = SaveStringFromText (gbp->taxName);
      
        /* make sure we use capitalization from list */
        if (gbp->orglist != NULL)
        {
          GetItemParams (gbp->orglist, 1, NULL, &num, NULL, NULL, NULL);
          if (num > 0) {
            row = FindTaxText (orp->taxname, num);
            if (row > 0 && row <= num) {
              CopyStrFromTaxPtr (txt, sizeof (txt) - 2, row, 1);
              if (StringICmp (txt, orp->taxname) == 0
                  && StringCmp (txt, orp->taxname) != 0) {
                orp->taxname = MemFree (orp->taxname);
                orp->taxname = StringSave (txt);
              }
            }
          }
        }
        

        /*
        orp->common = SaveStringFromText (gbp->commonName);
        */
        GetTitle (gbp->commonName, str, sizeof (str) - 1);
        TrimSpacesAroundString (str);
        if (! StringHasNoText (str)) {
          orp->common = StringSave (str);
        }
        orp->db = DialogToPointer (gbp->db);
        orp->syn = DialogToPointer (gbp->syn);
        orp->mod = DialogToPointer (gbp->mod);
        onp = OrgNameNew ();
        orp->orgname = onp;
        if (onp != NULL) {
          /* retrieve unaltered orgname data */
          onp->choice = gbp->orgname_choice;
          onp->data = gbp->orgname_data;
          gbp->orgname_choice = 0;
          gbp->orgname_data = NULL;

          if (gbp->simplecode != NULL) {
            vnp = DialogToPointer (gbp->genome);
            if (vnp != NULL) {
              genome = GenomeFromSrcLoc (vnp->choice);
              vnp = ValNodeFreeData (vnp);
              if (genome == 4 || genome == 5) {
                onp->mgcode = gcIndexToId [GetValue (gbp->simplecode)];
                onp->gcode = gcIndexToId [GetValue (gbp->gcode)];
              } else {
                onp->gcode = gcIndexToId [GetValue (gbp->simplecode)];
                onp->mgcode = gcIndexToId [GetValue (gbp->mgcode)];
              }
            }
          } else {
            onp->gcode = gcIndexToId [GetValue (gbp->gcode)];
            onp->mgcode = gcIndexToId [GetValue (gbp->mgcode)];
          }
          if (! TextHasNoText (gbp->gbDiv)) {
            onp->div = SaveStringFromText (gbp->gbDiv);
          }
          if (gbp->lineage == NULL && gbp->taxID != 0) {
            ProgramPath (str, sizeof (str));
            ptr = StringRChr (str, DIRDELIMCHR);
            if (ptr != NULL) {
              *ptr = '\0';
              FileBuildPath (str, NULL, "lineages.txt");
              f = FileOpen (str, "r");
              if (f == NULL) {
                if (GetAppParam ("NCBI", "NCBI", "DATA", "", str, sizeof (str))) {
                  FileBuildPath (str, NULL, "lineages.txt");
                  f = FileOpen (str, "r");
                }
              }
              if (f != NULL) {
                if (FileGets (str, sizeof (str), f) != NULL) {
                  goOn = (Boolean) (FileGets (str, sizeof (str), f) != NULL);
                  while (goOn) {
                    ptr = StringChr (str, '\t');
                    if (ptr != NULL) {
                      *ptr = '\0';
                      if (StrToLong (str, &taxID) && taxID == gbp->taxID) {
                        ptr++;
                        chptr = ptr;
                        ch = *chptr;
                        while (ch != '\0' && ch != '\r' && ch != '\n') {
                          chptr++;
                          ch = *chptr;
                        }
                        *chptr = '\0';
                        onp->lineage = StringSave (ptr);
                        goOn = FALSE;
                      }
                    }
                    goOn = (Boolean) (goOn && (FileGets (str, sizeof (str), f) != NULL));
                  }
                }
                FileClose (f);
              }
            }
          } else if (! TextHasNoText (gbp->lineage)) {
            onp->lineage = SaveStringFromTextAndStripNewlines (gbp->lineage);
          }
          onp->mod = DialogToPointer (gbp->orgmod_val_dlg);
          GetTitle (gbp->gbacr, buf, sizeof (buf) - 1);
          if (! StringHasNoText (buf)) {
            mod = OrgModNew ();
            if (onp->mod == NULL) {
              onp->mod = mod;
            } else {
              tmpmod = onp->mod;
              while (tmpmod->next != NULL) {
                tmpmod = tmpmod->next;
              }
              tmpmod->next = mod;
            }
            if (mod != NULL) {
              mod->subtype = 32;
              mod->subname = StringSave (buf);
            }
          }
          GetTitle (gbp->gbana, buf, sizeof (buf) - 1);
          if (! StringHasNoText (buf)) {
            mod = OrgModNew ();
            if (onp->mod == NULL) {
              onp->mod = mod;
            } else {
              tmpmod = onp->mod;
              while (tmpmod->next != NULL) {
                tmpmod = tmpmod->next;
              }
              tmpmod->next = mod;
            }
            if (mod != NULL) {
              mod->subtype = 33;
              mod->subname = StringSave (buf);
            }
          }
          GetTitle (gbp->gbsyn, buf, sizeof (buf) - 1);
          if (! StringHasNoText (buf)) {
            mod = OrgModNew ();
            if (onp->mod == NULL) {
              onp->mod = mod;
            } else {
              tmpmod = onp->mod;
              while (tmpmod->next != NULL) {
                tmpmod = tmpmod->next;
              }
              tmpmod->next = mod;
            }
            if (mod != NULL) {
              mod->subtype = 34;
              mod->subname = StringSave (buf);
            }
          }
          if (gbp->stripOldName && onp->mod != NULL) {
            prevmod = (OrgModPtr PNTR) &(onp->mod);
            tmpmod = onp->mod;
            while (tmpmod != NULL) {
              nextmod = tmpmod->next;
              if (tmpmod->subtype == 254) {
                *(prevmod) = tmpmod->next;
                tmpmod->next = NULL;
                OrgModFree (tmpmod);
              } else {
                prevmod = (OrgModPtr PNTR) &(tmpmod->next);
              }
              tmpmod = nextmod;
            }
          }
          if (onp->lineage == NULL && onp->mod == NULL &&
              onp->gcode == 0 && onp->mgcode == 0) {
            orp->orgname = OrgNameFree (orp->orgname);
          }
        }
      }

      biop->subtype = DialogToPointer (gbp->subsrc_val_dlg);

      RemoveTextFromTextFreeSubSourceModifiers (biop, NULL);     

      /* if we find plasmid-name on a location that cannot have
       * plasmids, change the location to plasmid */
      if (biop->genome != GENOME_mitochondrion
          && biop->genome != GENOME_chloroplast
          && biop->genome != GENOME_kinetoplast
          && biop->genome != GENOME_chromoplast
          && biop->genome != GENOME_plastid
          && biop->genome != GENOME_apicoplast
          && biop->genome != GENOME_leucoplast
          && biop->genome != GENOME_proplastid)
      {
        tmpssp = biop->subtype;
        while (tmpssp != NULL)
        {
          if (tmpssp->subtype == SUBSRC_plasmid_name) {
            biop->genome = GENOME_plasmid;
            break;
          }
          tmpssp = tmpssp->next;
        }
      }
 
      if (orp != NULL) {
        if (orp->taxname == NULL && orp->common == NULL &&
            orp->mod == NULL && orp->db == NULL && orp->syn == NULL &&
            orp->orgname == NULL && biop->subtype == NULL) {
          biop = BioSourceFree (biop);
        }
      } else {
        biop = BioSourceFree (biop);
      }
    }
  }
  return (Pointer) biop;
}

static ValNodePtr TestGenBioDialog (DialoG d)

{
  Char           comm [64];
  GenBioPagePtr  gbp;
  ValNodePtr     head;

  head = NULL;
  gbp = (GenBioPagePtr) GetObjectExtra (d);
  if (gbp != NULL) {
    GetTitle (gbp->commonName, comm, sizeof (comm) - 1);
    if (TextHasNoText (gbp->taxName) && StringHasNoText (comm) ) {
      head = AddStringToValNodeChain (head, "organism name", 1);
    }
  }
  return head;
}

static void BioSourceMessage (DialoG d, Int2 mssg)

{
  GenBioPagePtr  gbp;

  gbp = (GenBioPagePtr) GetObjectExtra (d);
  if (gbp != NULL) {
    switch (mssg) {
      case VIB_MSG_INIT :
        UpdateDocument (gbp->orglist, 0, 0);
        AdjustDocScroll (gbp->orglist);
        break;
      case VIB_MSG_ENTER :
        Select (gbp->taxName);
        break;
      default :
        break;
    }
  }
}


extern PopuP ReplaceBioSourceGencodePopup (DialoG d, PopuP gencode);
PopuP ReplaceBioSourceGencodePopup (DialoG d, PopuP gencode)

{
  GenBioPagePtr  gbp;
  PopuP          orig_gencode = NULL;

  gbp = (GenBioPagePtr) GetObjectExtra (d);
  if (gbp != NULL) {
    orig_gencode = gbp->simplecode;
    gbp->simplecode = gencode;
  }
  return orig_gencode;
}


NLM_EXTERN ValNodePtr GetLocListForBioSource (BioSourcePtr biop)
{
  ValNodePtr loc_list, vnp, vnp_prev = NULL, vnp_next;
  Boolean    indexerVersion;

  indexerVersion = (Boolean) (GetAppProperty ("InternalNcbiSequin") != NULL);

  loc_list = GetLocationList (FALSE);

  /* remove chromosome, transposon, insertion seq if not already on biosource */
  
  vnp = loc_list;
  while (vnp != NULL) {
    vnp_next = vnp->next;
    if ((vnp->choice == Source_location_chromosome && !indexerVersion
          && (biop == NULL || biop->genome != GENOME_chromosome))
        || (vnp->choice == Source_location_virion
            && (biop == NULL || biop->genome != GENOME_virion))
        || (vnp->choice == Source_location_transposon
            && (biop == NULL || biop->genome != GENOME_transposon))
        || (vnp->choice == Source_location_insertion_seq
            && (biop == NULL || biop->genome != GENOME_insertion_seq))) {
      if (vnp_prev == NULL) { 
        loc_list = vnp->next;
      } else {
        vnp_prev->next= vnp->next;
      }
      vnp->next = NULL;
      vnp = ValNodeFreeData (vnp);
    } else {
      vnp_prev = vnp;
    }
    vnp = vnp_next;
  }
  return loc_list;
}


extern DialoG CreateSimpleBioSourceDialog (GrouP h, CharPtr title)

{
  GrouP          f;
  GrouP          g;
  GenBioPagePtr  gbp;
  Int2           height;
  GrouP          m;
  GrouP          p;
  GrouP          q;
  RecT           r;
  GrouP          s;
  GrouP          x;
  ValNode        vn;

  p = HiddenGroup (h, 1, 0, NULL);
  SetGroupSpacing (p, 10, 10);

  gbp = (GenBioPagePtr) MemNew (sizeof (GenBioPage));
  if (gbp != NULL) {

    SetObjectExtra (p, gbp, StdCleanupExtraProc);
    gbp->dialog = (DialoG) p;
    gbp->todialog = BioSourcePtrToGenBioPage;
    gbp->fromdialog = GenBioPageToBioSourcePtr;
    gbp->dialogmessage = BioSourceMessage;
    gbp->testdialog = TestGenBioDialog;

    if (title != NULL && title [0] != '\0') {
      s = NormalGroup (p, 0, -2, title, systemFont, NULL);
    } else {
      s = HiddenGroup (p, 0, -2, NULL);
    }
    SetGroupSpacing (s, 10, 10);

    m = HiddenGroup (s, 0, 0, NULL);

    q = HiddenGroup (m, -1, 0, NULL);
    SetGroupSpacing (q, 10, 10);

    gbp->origTaxName = NULL;
    gbp->stripOldName = FALSE;

    g = HiddenGroup (q, 1, 0, NULL);
    /*
    StaticPrompt (g, "Organism", 0, 0, programFont, 'c');
    */
    f = HiddenGroup (g, 2, 0, NULL);
    StaticPrompt (f, "Scientific Name", 0, dialogTextHeight, programFont, 'l');
    gbp->taxName = DialogText (f, "", 20, TaxNameText);
    SetObjectExtra (gbp->taxName, gbp, NULL);
    StaticPrompt (f, "Common Name", 0, dialogTextHeight, programFont, 'l');
    /*
    gbp->commonName = DialogText (f, "", 10, CommonNameText);
    */
    gbp->commonName = (Handle) StaticPrompt (f, "", 10 * stdCharWidth, dialogTextHeight, systemFont, 'l');
    SetObjectExtra (gbp->commonName, gbp, NULL);
    StaticPrompt (f, "", 0, dialogTextHeight, programFont, 'l');
    f = HiddenGroup (g, 1, 0, NULL);
    SelectFont (programFont);
    height = LineHeight ();
    SelectFont (systemFont);
    gbp->orglist = DocumentPanel (f, stdCharWidth * 25, height * 6);
    SetObjectExtra (gbp->orglist, gbp, NULL);
    SetDocAutoAdjust (gbp->orglist, TRUE);
    orgListCol [0].pixWidth = screenRect.right - screenRect.left;
    AppendItem (gbp->orglist, AllButFirstLinePrtProc, orgTxtPtr, FALSE, orgNum,
                &orgListPar, orgListCol, programFont);
    SetDocAutoAdjust (gbp->orglist, TRUE);
    SetDocProcs (gbp->orglist, ClickTaxName, NULL, ReleaseTaxName, NULL);
    SetDocShade (gbp->orglist, NULL, NULL, HighlightTaxName, NULL);

    AlignObjects (ALIGN_RIGHT, (HANDLE) gbp->taxName, (HANDLE) gbp->commonName,
                  (HANDLE) gbp->orglist, NULL);

    g = HiddenGroup (q, -1, 0, NULL);
    f = HiddenGroup (g, 3, 0, NULL);
    StaticPrompt (f, "Location of Sequence",
                  0, popupMenuHeight, programFont, 'l');

    gbp->genome = ValNodeSelectionDialogExEx (f, GetLocListForBioSource (NULL), SHORT_SELECTION_LIST, ValNodeStringName,
                                           ValNodeSimpleDataFree, ValNodeStringCopy,
                                           ValNodeChoiceMatch, "location", 
                                           NULL, NULL, FALSE, FALSE, TRUE, NULL);
    vn.choice = Source_location_genomic;
    vn.data.ptrvalue = NULL;
    vn.next = NULL;
    PointerToDialog (gbp->genome, &vn);
    ObjectRect (gbp->orglist, &r);
    MultiLinePrompt (g, useGenomicText, r.right - r.left - 2, programFont);

    x = HiddenGroup (q, 0, 0, NULL);

    f = HiddenGroup (x, 2, 0, NULL);
    StaticPrompt (f, "Genetic Code for Translation", 0, popupMenuHeight, programFont, 'l');
    gbp->simplecode = PopupList (f, TRUE, NULL);
    PopulateGeneticCodePopup (gbp->simplecode);
    SetValue (gbp->simplecode, 1);
    gbp->gbDiv = DialogText (x, "", 4, NULL);
    Hide (gbp->gbDiv);

/* superimpose two hidden genetic code controls to save both in resulting biosource */

    gbp->gcode = PopupList (x, TRUE, NULL);
    PopulateGeneticCodePopup (gbp->gcode);
    SetValue (gbp->gcode, 1);
    gbp->mgcode = PopupList (x, TRUE, NULL);
    PopulateGeneticCodePopup (gbp->mgcode);
    SetValue (gbp->mgcode, 1);
    Hide (gbp->gcode);
    Hide (gbp->mgcode);

    SelectFont (systemFont);
  }

  return (DialoG) p;
}

static void OrgModPtrToOrgmodDialog (DialoG d, Pointer data)

{
  ValNodePtr  head;
  Int2        j;
  size_t      len;
  OrgModPtr   list;
  CharPtr     str;
  TagListPtr  tlp;
  Char        tmp [16];
  ValNodePtr  vnp;

  tlp = (TagListPtr) GetObjectExtra (d);
  list = (OrgModPtr) data;
  if (tlp != NULL) {
    head = NULL;
    while (list != NULL) {
      if (list->subname != NULL && list->subtype != 255 &&
          list->subtype != 32 && list->subtype != 33 && list->subtype != 34) {
        vnp = ValNodeNew (head);
        if (head == NULL) {
          head = vnp;
        }
        if (vnp != NULL) {
          sprintf (tmp, "%d", (int) list->subtype);
          len = StringLen (tmp) + StringLen (list->subname);
          str = MemNew (len + 4);
          if (str != NULL) {
            StringCpy (str, tmp);
            StringCat (str, "\t");
            StringCat (str, list->subname);
            StringCat (str, "\n");
          }
          vnp->data.ptrvalue = str;
        }
      }
      list = list->next;
    }
    SendMessageToDialog (tlp->dialog, VIB_MSG_RESET);
    tlp->vnp = head;
    SendMessageToDialog (tlp->dialog, VIB_MSG_REDRAW);
    for (j = 0, vnp = tlp->vnp; vnp != NULL; j++, vnp = vnp->next) {
    }
    tlp->max = MAX ((Int2) 0, (Int2) (j - tlp->rows + 1));
    CorrectBarMax (tlp->bar, tlp->max);
    CorrectBarPage (tlp->bar, tlp->rows - 1, tlp->rows - 1);
  }
}

static Pointer OrgmodDialogToOrgModPtr (DialoG d)

{
  Char        ch;
  OrgModPtr   head;
  Int2        j;
  Int2        len;
  OrgModPtr   omp;
  OrgModPtr   omplast;
  Boolean     okay;
  CharPtr     str;
  TagListPtr  tlp;
  CharPtr     tmp;
  int         val;
  ValNodePtr  vnp;

  head = NULL;
  tlp = (TagListPtr) GetObjectExtra (d);
  if (tlp != NULL && tlp->vnp != NULL) {
    omp = NULL;
    omplast = NULL;
    for (vnp = tlp->vnp; vnp != NULL; vnp = vnp->next) {
      str = (CharPtr) vnp->data.ptrvalue;
      okay = FALSE;
      len = StringLen (str);
      for (j = 0; j < len; j++) {
        ch = str [j];
        if (ch != ' ' && ch != '\t' && ch != '\n') {
          okay = TRUE;
        }
      }
      if (okay) {
        tmp = ExtractTagListColumn ((CharPtr) vnp->data.ptrvalue, 0);
        if (tmp != NULL && sscanf (tmp, "%d", &val) == 1 && val != 0) {
          MemFree (tmp);
          tmp = ExtractTagListColumn ((CharPtr) vnp->data.ptrvalue, 1);
          if (! StringHasNoText (tmp)) {
            omp = OrgModNew ();
            if (omplast == NULL) {
              head = omp;
            } else {
              omplast->next = omp;
            }
            omplast = omp;
            if (omp != NULL) {
              omp->subtype = (Uint1) val;
              omp->subname = tmp;
            }
          } else {
            MemFree (tmp);
          }
        } else {
          MemFree (tmp);
        }
      }
    }
  }
  return (Pointer) head;
}

static void SubSourcePtrToSubsourceDialog (DialoG d, Pointer data)

{
  ValNodePtr    head;
  Int2          j;
  size_t        len;
  SubSourcePtr  list;
  CharPtr       str;
  TagListPtr    tlp;
  Char          tmp [16];
  ValNodePtr    vnp;

  tlp = (TagListPtr) GetObjectExtra (d);
  list = (SubSourcePtr) data;
  if (tlp != NULL) {
    head = NULL;
    while (list != NULL) {
      if (list->name != NULL && list->subtype != 255) {
        vnp = ValNodeNew (head);
        if (head == NULL) {
          head = vnp;
        }
        if (vnp != NULL) {
          sprintf (tmp, "%d", (int) list->subtype);
          len = StringLen (tmp) + StringLen (list->name);
          str = MemNew (len + 4);
          if (str != NULL) {
            StringCpy (str, tmp);
            StringCat (str, "\t");
            StringCat (str, list->name);
            StringCat (str, "\n");
          }
          vnp->data.ptrvalue = str;
        }
      }
      list = list->next;
    }
    SendMessageToDialog (tlp->dialog, VIB_MSG_RESET);
    tlp->vnp = head;
    SendMessageToDialog (tlp->dialog, VIB_MSG_REDRAW);
    for (j = 0, vnp = tlp->vnp; vnp != NULL; j++, vnp = vnp->next) {
    }
    tlp->max = MAX ((Int2) 0, (Int2) (j - tlp->rows + 1));
    CorrectBarMax (tlp->bar, tlp->max);
    CorrectBarPage (tlp->bar, tlp->rows - 1, tlp->rows - 1);
  }
}

typedef struct fixmodifiertextform
{
  WindoW     w;
  Boolean    done;
  Boolean    move_to_text;
  Boolean    remove;
  Boolean    do_all;
} FixModifierTextFormData, PNTR FixModifierTextFormPtr;

static void FixModifierTextMove (ButtoN b)
{
  FixModifierTextFormPtr fp;
  
  fp = (FixModifierTextFormPtr) GetObjectExtra (b);
  if (fp == NULL) return;
  
  Remove (fp->w);
  fp->remove = FALSE;
  fp->move_to_text = TRUE;
  fp->do_all = FALSE;
  fp->done = TRUE;    
}

static void FixModifierTextMoveAll (ButtoN b)
{
  FixModifierTextFormPtr fp;
  
  fp = (FixModifierTextFormPtr) GetObjectExtra (b);
  if (fp == NULL) return;
  
  Remove (fp->w);
  fp->remove = FALSE;
  fp->move_to_text = TRUE;
  fp->do_all = TRUE;
  fp->done = TRUE;    
}

static void FixModifierTextRemove (ButtoN b)
{
  FixModifierTextFormPtr fp;
  
  fp = (FixModifierTextFormPtr) GetObjectExtra (b);
  if (fp == NULL) return;
  
  Remove (fp->w);
  fp->remove = TRUE;
  fp->move_to_text = FALSE;
  fp->do_all = FALSE;
  fp->done = TRUE;    
}

static void FixModifierTextRemoveAll (ButtoN b)
{
  FixModifierTextFormPtr fp;
  
  fp = (FixModifierTextFormPtr) GetObjectExtra (b);
  if (fp == NULL) return;
  
  Remove (fp->w);
  fp->remove = TRUE;
  fp->move_to_text = FALSE;
  fp->do_all = TRUE;
  fp->done = TRUE;    
}

extern ModTextFixPtr ModTextFixNew (void)
{
  ModTextFixPtr tfp;
  
  tfp = (ModTextFixPtr) MemNew (sizeof (ModTextFixData));
  if (tfp == NULL) return NULL;
  tfp->remove_this = FALSE;
  tfp->move_this = FALSE;
  tfp->remove_all_germline = FALSE;
  tfp->remove_all_transgenic = FALSE;
  tfp->remove_all_environmental = FALSE;
  tfp->remove_all_rearranged = FALSE;
  tfp->remove_all_metagenomic = FALSE;
  tfp->move_all_germline = FALSE;
  tfp->move_all_transgenic = FALSE;
  tfp->move_all_environmental = FALSE;
  tfp->move_all_rearranged = FALSE;
  tfp->move_all_metagenomic = FALSE;
  return tfp;
}

static void 
GetModifierTextFix (ModTextFixPtr tfp, Uint1 subtype, CharPtr txt)
{
  GrouP  g, c, t;
  ButtoN b;
  FixModifierTextFormData fd;
  CharPtr prompt_fmt = "You have text (%s) in %s modifier field.";
  CharPtr prompt_str = NULL;
  
  if (tfp == NULL) return;
  switch (subtype)
  {
      case SUBSRC_rearranged:
        if (tfp->remove_all_rearranged)
        {
            tfp->remove_this = TRUE;
            tfp->move_this = FALSE;
            return;
        }
        else if (tfp->move_all_rearranged)
        {
            tfp->move_this = TRUE;
            tfp->remove_this = FALSE;
            return;
        }
        break;
      case SUBSRC_transgenic:
        if (tfp->remove_all_transgenic)
        {
            tfp->remove_this = TRUE;
            tfp->move_this = FALSE;
            return;
        }
        else if (tfp->move_all_transgenic)
        {
            tfp->move_this = TRUE;
            tfp->remove_this = FALSE;
            return;
        }
        break;
      case SUBSRC_germline:
        if (tfp->remove_all_germline)
        {
            tfp->remove_this = TRUE;
            tfp->move_this = FALSE;
            return;
        }
        else if (tfp->move_all_germline)
        {
            tfp->move_this = TRUE;
            tfp->remove_this = FALSE;
            return;
        }
        break;
      case SUBSRC_environmental_sample:
        if (tfp->remove_all_environmental)
        {
            tfp->remove_this = TRUE;
            tfp->move_this = FALSE;
            return;
        }
        else if (tfp->move_all_environmental)
        {
            tfp->move_this = TRUE;
            tfp->remove_this = FALSE;
            return;
        }
        break;
      case SUBSRC_metagenomic:
        if (tfp->remove_all_metagenomic)
        {
            tfp->remove_this = TRUE;
            tfp->move_this = FALSE;
            return;
        }
        else if (tfp->move_all_metagenomic)
        {
            tfp->move_this = TRUE;
            tfp->remove_this = FALSE;
            return;
        }
        break;
      default:
        break;
  }

  fd.w = ModalWindow(-20, -13, -10, -10, NULL);
  g = HiddenGroup(fd.w, -1, 0, NULL);
  
  prompt_str = (CharPtr) MemNew (sizeof (Char) * (StringLen (prompt_fmt) + StringLen (txt)
                                  + StringLen ("an environmental sample")));
  if (prompt_str == NULL) return;
  switch (subtype)
  {
      case SUBSRC_rearranged:
        sprintf (prompt_str, prompt_fmt, txt, "a rearranged");
        break;
      case SUBSRC_germline:
        sprintf (prompt_str, prompt_fmt, txt, "a germline");
        break;
      case SUBSRC_transgenic:
        sprintf (prompt_str, prompt_fmt, txt, "a transgenic");
        break;
      case SUBSRC_environmental_sample:
        sprintf (prompt_str, prompt_fmt, txt, "an environmental sample");
        break;
      case SUBSRC_metagenomic:
        sprintf (prompt_str, prompt_fmt, txt, "a metagenomic");
        break;
  }
  
  t = HiddenGroup (g, 1, 0, NULL);
  StaticPrompt (t, prompt_str, 0, dialogTextHeight, programFont, 'l');
  StaticPrompt (t, "This text will never be displayed in your GenBank record.", 0, dialogTextHeight, programFont, 'l');
  StaticPrompt (t, "Do you want to move this text to a note or remove it?", 0, dialogTextHeight, programFont, 'l');
  
  c = HiddenGroup (g, 4, 0, NULL);
  b = PushButton (c, "Move to note", FixModifierTextMove);
  SetObjectExtra (b, &fd, NULL);  
  b = PushButton (c, "Move all to note", FixModifierTextMoveAll);
  SetObjectExtra (b, &fd, NULL);
  b = PushButton (c, "Remove", FixModifierTextRemove);
  SetObjectExtra (b, &fd, NULL);
  b = PushButton (c, "Remove all", FixModifierTextRemoveAll);
  SetObjectExtra (b, &fd, NULL);
  AlignObjects (ALIGN_CENTER, (HANDLE) t, (HANDLE) c, NULL);
  
  Show(fd.w); 
  Select (fd.w);
  fd.done = FALSE;
  while (!fd.done)
  {
    ProcessExternalEvent ();
    Update ();
  }
  ProcessAnEvent ();

  if (fd.remove)
  {
      tfp->remove_this = TRUE;
      if (fd.do_all)
      {
        switch (subtype)
        {
            case SUBSRC_rearranged:
              tfp->remove_all_rearranged = TRUE;
              tfp->move_all_rearranged = FALSE;
              break;
            case SUBSRC_transgenic:
              tfp->remove_all_transgenic = TRUE;
              tfp->move_all_transgenic = FALSE;
              break;
            case SUBSRC_germline:
              tfp->remove_all_germline = TRUE;
              tfp->move_all_germline = FALSE;
              break;
            case SUBSRC_environmental_sample:
              tfp->remove_all_environmental = TRUE;
              tfp->move_all_environmental = FALSE;
              break;
            case SUBSRC_metagenomic:
              tfp->remove_all_metagenomic = TRUE;
              tfp->move_all_metagenomic = FALSE;
              break;
        }
      }
  }
  else if (fd.move_to_text)
  {
      tfp->move_this = TRUE;
      if (fd.do_all)
      {
        switch (subtype)
        {
            case SUBSRC_rearranged:
              tfp->remove_all_rearranged = FALSE;
              tfp->move_all_rearranged = TRUE;
              break;
            case SUBSRC_transgenic:
              tfp->remove_all_transgenic = FALSE;
              tfp->move_all_transgenic = TRUE;
              break;
            case SUBSRC_germline:
              tfp->remove_all_germline = FALSE;
              tfp->move_all_germline = TRUE;
              break;
            case SUBSRC_environmental_sample:
              tfp->remove_all_environmental = FALSE;
              tfp->move_all_environmental = TRUE;
              break;
            case SUBSRC_metagenomic:
              tfp->remove_all_metagenomic = FALSE;
              tfp->move_all_metagenomic = TRUE;
              break;
        }
      }
  }
}

extern void RemoveTextFromTextFreeSubSourceModifiers (BioSourcePtr biop, Pointer userdata)
{
  SubSourcePtr ssp;
  SubSourcePtr note_ssp = NULL;
  Int4         len;
  CharPtr      new_note;
  ModTextFixPtr   tfp;
  
  if (biop == NULL || biop->subtype == NULL) return;
  
  if (userdata == NULL)
  {
    tfp = ModTextFixNew();
    if (tfp == NULL) return;
  }
  else
  {
      tfp = (ModTextFixPtr) userdata;
  }
  
  for (ssp = biop->subtype; ssp != NULL; ssp = ssp->next)
  {
      tfp->move_this = FALSE;
      tfp->remove_this = FALSE;
      if ((ssp->subtype == SUBSRC_germline
          || ssp->subtype == SUBSRC_transgenic
          || ssp->subtype == SUBSRC_rearranged
          || ssp->subtype == SUBSRC_environmental_sample
          || ssp->subtype == SUBSRC_metagenomic)
          && ! StringHasNoText (ssp->name))
      {
       GetModifierTextFix (tfp, ssp->subtype, ssp->name);
        if (tfp->move_this)
        {
          /* if a note modifier is found, add this text to it, otherwise create a new
           * note modifier to hold this text.
           */
          if (note_ssp == NULL)
          {
            for (note_ssp = biop->subtype; note_ssp != NULL && note_ssp->subtype != 255; note_ssp = note_ssp->next)
            {    
            }
          }
          if (note_ssp == NULL)
          {
              note_ssp = SubSourceNew ();
              if (note_ssp != NULL)
              {
                note_ssp->subtype = 255;
                note_ssp->name = ssp->name;
                ssp->name = StringSave ("");
                note_ssp->next = ssp->next;
                ssp->next = note_ssp;
              }
          }
          else if (StringHasNoText (note_ssp->name))
          {
          note_ssp->name = MemFree (note_ssp->name);
          note_ssp->name = ssp->name;
          ssp->name = StringSave ("");                
          }
          else
          {
              len = StringLen (note_ssp->name) + StringLen (ssp->name) + 3;
              new_note = (CharPtr) MemNew (len * sizeof (Char));
              if (new_note != NULL)
              {
                StringCpy (new_note, note_ssp->name);
                StringCat (new_note, "; ");
                StringCat (new_note, ssp->name);
                note_ssp->name = MemFree (note_ssp->name);
                note_ssp->name = new_note;
                ssp->name = MemFree (ssp->name);
                ssp->name = StringSave ("");
              }
          }
        }
        else if (tfp->remove_this)
        {
          ssp->name = MemFree (ssp->name);
          ssp->name = StringSave ("");
        }
      }
  }
  if (userdata == NULL)
  {
      MemFree (tfp);
  }
}

static Pointer SubsourceDialogToSubSourcePtr (DialoG d)

{
  Char          ch;
  SubSourcePtr  head;
  Int2          j;
  Int2          len;
  SubSourcePtr  ssp;
  SubSourcePtr  ssplast;
  Boolean       okay;
  CharPtr       str;
  TagListPtr    tlp;
  CharPtr       tmp;
  int           val;
  ValNodePtr    vnp;

  head = NULL;
  tlp = (TagListPtr) GetObjectExtra (d);
  if (tlp != NULL && tlp->vnp != NULL) {
    ssp = NULL;
    ssplast = NULL;
    for (vnp = tlp->vnp; vnp != NULL; vnp = vnp->next) {
      str = (CharPtr) vnp->data.ptrvalue;
      okay = FALSE;
      len = StringLen (str);
      for (j = 0; j < len; j++) {
        ch = str [j];
        if (ch != ' ' && ch != '\t' && ch != '\n') {
          okay = TRUE;
        }
      }
      if (okay) {
        tmp = ExtractTagListColumn ((CharPtr) vnp->data.ptrvalue, 0);
        if (tmp != NULL && sscanf (tmp, "%d", &val) == 1 && val != 0) {
          MemFree (tmp);
          tmp = ExtractTagListColumn ((CharPtr) vnp->data.ptrvalue, 1);
          if ((val == SUBSRC_germline ||
               val == SUBSRC_rearranged ||
               val == SUBSRC_transgenic ||
               val == SUBSRC_environmental_sample ||
               val == SUBSRC_metagenomic) &&
              StringHasNoText (tmp)) {
            MemFree (tmp);
            tmp = StringSave ("");
          }
          if ((! StringHasNoText (tmp)) ||
              val == SUBSRC_germline ||
              val == SUBSRC_rearranged ||
              val == SUBSRC_transgenic ||
              val == SUBSRC_environmental_sample ||
              val == SUBSRC_metagenomic) {
            ssp = SubSourceNew ();
            if (ssplast == NULL) {
              head = ssp;
            } else {
              ssplast->next = ssp;
            }
            ssplast = ssp;
            if (ssp != NULL) {
              ssp->subtype = (Uint1) val;
              ssp->name = tmp;
              tmp = NULL;
            }
          }
          MemFree (tmp);
        }
      }
    }
  }
  return (Pointer) head;
}

Uint2 orgmod_widths [] = {
  0, 25
};

Uint2 subsource_widths [] = {
  0, 25
};

Uint2 orgmod_types [] = {
  TAGLIST_POPUP, TAGLIST_TEXT
};

Uint2 subsource_types [] = {
  TAGLIST_POPUP, TAGLIST_TEXT
};

/*
static CharPtr orgmod_extra_prompts [] = {
  "Additional", "Organism", "Information", NULL
};

static CharPtr subsource_extra_prompts [] = {
  "Additional", "Source", "Information", NULL
};
*/

static CharPtr orgTabs [] = {
  "Names", "Location", "Genetic Codes", "Lineage", NULL
};

static CharPtr modTabs [] = {
  "Source", "Organism", "GenBank", NULL
};

static CharPtr modTabsUns [] = {
  "Source", "Organism", "GenBank", "Unstructured", NULL
};

static CharPtr miscTabs1 [] = {
  "Cross-Refs", NULL
};

static CharPtr miscTabs2 [] = {
  "Synonyms", "Cross-Refs", NULL
};

static void ChangeOrgSubPage (VoidPtr data, Int2 newval, Int2 oldval)

{
  GenBioPagePtr  gbp;

  gbp = (GenBioPagePtr) data;
  if (gbp != NULL) {
    if (oldval >= 0 && oldval <= 3) {
      SafeHide (gbp->orgGrp [oldval]);
    }
    if (newval >= 0 && newval <= 3) {
      SafeShow (gbp->orgGrp [newval]);
    }
    Update ();
  }
}

static void ChangeModSubPage (VoidPtr data, Int2 newval, Int2 oldval)

{
  GenBioPagePtr  gbp;

  gbp = (GenBioPagePtr) data;
  if (gbp != NULL) {
    if (oldval >= 0 && oldval <= 3) {
      SafeHide (gbp->modGrp [oldval]);
    }
    if (newval >= 0 && newval <= 3) {
      SafeShow (gbp->modGrp [newval]);
    }
    Update ();
  }
}

static void ChangeMiscSubPage (VoidPtr data, Int2 newval, Int2 oldval)

{
  GenBioPagePtr  gbp;

  gbp = (GenBioPagePtr) data;
  if (gbp != NULL) {
    if (oldval >= 0 && oldval <= 1) {
      SafeHide (gbp->miscGrp [oldval]);
    }
    if (newval >= 0 && newval <= 1) {
      SafeShow (gbp->miscGrp [newval]);
    }
    Update ();
  }
}

static void LookupTheTaxonomyProc (ButtoN b)

{
  GenBioFormPtr  gfp;

  gfp = (GenBioFormPtr) GetObjectExtra (b);
  if (gfp != NULL && gfp->lookupTaxonomy != NULL) {
    if (gfp->lookupTaxonomy (gfp->input_entityID)) {
      Remove (gfp->form);
      ObjMgrSendMsg (OM_MSG_UPDATE, gfp->input_entityID, 0, 0);
      Update ();
    }
  }
}

static EnumFieldAssocPtr EnumListFromQualNameAssoc (Nlm_QualNameAssocPtr qp)
{
  EnumFieldAssocPtr eap;
  /* start num_qual at one to count terminator */
  Int4              i, num_qual = 1;

  for (i = 0; qp[i].name != NULL; i++) {
    num_qual++;
  }

  eap = (EnumFieldAssocPtr) MemNew (sizeof (EnumFieldAssoc) * num_qual);
  for (i = 0; qp[i].name != NULL; i++) {
    eap[i].name = qp[i].name;
    eap[i].value = qp[i].value;
  }
  eap[i].name = NULL;
  eap[i].value = 0;
  return eap;
}

static void CleanupBioSourceDialog (GraphiC g, VoidPtr data)

{
  GenBioPagePtr  gbp;

  gbp = (GenBioPagePtr) data;
  if (gbp != NULL) {
    gbp->origTaxName = MemFree (gbp->origTaxName);
    gbp->orgmod_alists[0] = MemFree (gbp->orgmod_alists[0]);
    gbp->subsource_alists[0] = MemFree (gbp->subsource_alists[0]);
    FreeGenBioOrgNameData (gbp);
  }
  StdCleanupExtraProc (g, data);
}


static DialoG CreateBioSourceDialog (GrouP h, CharPtr title, GrouP PNTR pages,
                                     BioSourcePtr biop, GenBioFormPtr gfp,
                                     Boolean diableFocusControl)

{
  ButtoN         b;
  GrouP          c;
  GrouP          f, f1, f2, f3;
  GrouP          g;
  GenBioPagePtr  gbp;
  Boolean        hasSynonyms;
  Int2           height;
  Char           just;
  GrouP          k;
  GrouP          m;
  OrgRefPtr      orp;
  GrouP          p;
  PrompT         ppt;
  GrouP          q;
  RecT           r;
  GrouP          s;
  Boolean        showUnstructMods;
  GrouP          t;
  CharPtr PNTR   tabs = NULL;
  DialoG         tbs;
  PrompT         y;
  Int2           z;
  Boolean        indexerVersion;
  Boolean        has_discontinued, has_discouraged;
  ValNode        vn;

  p = HiddenGroup (h, 1, 0, NULL);
  SetGroupSpacing (p, 10, 10);

  gbp = (GenBioPagePtr) MemNew (sizeof (GenBioPage));
  if (gbp != NULL && pages != NULL) {

    SetObjectExtra (p, gbp, CleanupBioSourceDialog);
    gbp->dialog = (DialoG) p;
    gbp->todialog = BioSourcePtrToGenBioPage;
    gbp->fromdialog = GenBioPageToBioSourcePtr;
    gbp->dialogmessage = BioSourceMessage;
    gbp->testdialog = TestGenBioDialog;
    gbp->orgname_choice = 0;
    gbp->orgname_data = NULL;

    /* set up subsource and orgmod alists */
    indexerVersion = (Boolean) (GetAppProperty ("InternalNcbiSequin") != NULL);
    BioSourceHasOldOrgModQualifiers (biop, &has_discouraged, &has_discontinued);
    gbp->orgmod_alists[0] = GetModifiersEnum (FALSE, TRUE, has_discouraged || indexerVersion, has_discontinued);
    gbp->orgmod_alists[1] = NULL;
    BioSourceHasOldSubSourceQualifiers (biop, &has_discouraged, &has_discontinued);
    gbp->subsource_alists[0] = GetModifiersEnum (TRUE, FALSE, has_discouraged || indexerVersion, has_discontinued);
    gbp->subsource_alists[1] = NULL;

    if (title != NULL && title [0] != '\0') {
      s = NormalGroup (p, 0, -2, title, systemFont, NULL);
    } else {
      s = HiddenGroup (p, 0, -2, NULL);
    }
    SetGroupSpacing (s, 10, 10);

    m = HiddenGroup (s, 0, 0, NULL);

    pages [0] = HiddenGroup (m, -1, 0, NULL);
    SetGroupSpacing (pages [0], 10, 10);

    tbs = CreateFolderTabs (pages [0], orgTabs, 0, 0, 0,
                            PROGRAM_FOLDER_TAB,
                            ChangeOrgSubPage, (Pointer) gbp);
    k = HiddenGroup (pages [0], 0, 0, NULL);

    gbp->orgGrp [0] = HiddenGroup (k, -1, 0, NULL);
    SetGroupSpacing (gbp->orgGrp [0], 10, 10);

    gbp->origTaxName = NULL;
    gbp->stripOldName = FALSE;

    g = HiddenGroup (gbp->orgGrp [0], 1, 0, NULL);
    SetGroupSpacing (g, 10, 10);
    f = HiddenGroup (g, 2, 0, NULL);
    SetGroupSpacing (f, 3, 5);
    StaticPrompt (f, "Scientific Name", 0, dialogTextHeight, programFont, 'l');
    gbp->taxName = DialogText (f, "", 20, TaxNameText);
    SetObjectExtra (gbp->taxName, gbp, NULL);
    StaticPrompt (f, "Common Name", 0, dialogTextHeight, programFont, 'l');
    gbp->commonName = (Handle) DialogText (f, "", 10, CommonNameText);
    /*
    gbp->commonName = StaticPrompt (f, "", 10 * stdCharWidth, dialogTextHeight, systemFont, 'l');
    */
    SetObjectExtra (gbp->commonName, gbp, NULL);
    StaticPrompt (f, "", 0, dialogTextHeight, programFont, 'l');
    f = HiddenGroup (g, 1, 0, NULL);
    SelectFont (programFont);
    height = LineHeight ();
    SelectFont (systemFont);
    gbp->orglist = DocumentPanel (f, stdCharWidth * 25, height * 6);
    SetObjectExtra (gbp->orglist, gbp, NULL);
    SetDocAutoAdjust (gbp->orglist, FALSE);
    orgListCol [0].pixWidth = screenRect.right - screenRect.left;
    AppendItem (gbp->orglist, AllButFirstLinePrtProc, orgTxtPtr, FALSE, orgNum,
                &orgListPar, orgListCol, programFont);
    SetDocAutoAdjust (gbp->orglist, TRUE);
    SetDocProcs (gbp->orglist, ClickTaxName, NULL, ReleaseTaxName, NULL);
    SetDocShade (gbp->orglist, NULL, NULL, HighlightTaxName, NULL);

    AlignObjects (ALIGN_RIGHT, (HANDLE) gbp->taxName, (HANDLE) gbp->commonName,
                  (HANDLE) gbp->orglist, NULL);

    gbp->orgGrp [1] = HiddenGroup (k, -1, 0, NULL);
    SetGroupSpacing (gbp->orgGrp [1], 10, 10);

    g = HiddenGroup (gbp->orgGrp [1], -1, 0, NULL);
    SelectFont (programFont);
    f = HiddenGroup (g, 3, 0, NULL);
    StaticPrompt (f, "Location of Sequence",
                  0, popupMenuHeight, programFont, 'l');

    gbp->genome = ValNodeSelectionDialogExEx (f, GetLocListForBioSource (biop), SHORT_SELECTION_LIST, ValNodeStringName,
                                           ValNodeSimpleDataFree, ValNodeStringCopy,
                                           ValNodeChoiceMatch, "location", 
                                           NULL, NULL, FALSE, FALSE, TRUE, NULL);
    vn.choice = Source_location_genomic;
    vn.data.ptrvalue = NULL;
    vn.next = NULL;
    PointerToDialog (gbp->genome, &vn);

    ObjectRect (gbp->orglist, &r);
    MultiLinePrompt (g, useGenomicText, r.right - r.left - 2, programFont);

    f = HiddenGroup (gbp->orgGrp [1], 3, 0, NULL);
    StaticPrompt (f, "Origin of Sequence",
                  0, popupMenuHeight, programFont, 'l');
    gbp->origin = PopupList (f, TRUE, NULL);
    SetObjectExtra (gbp->origin, gbp, NULL);
    InitEnumPopup (gbp->origin, biosource_origin_alist, NULL);
    SetValue (gbp->origin, 0);

    gbp->is_focus = CheckBox (gbp->orgGrp [1], "Biological focus (if multiple source features)", NULL);
    if (diableFocusControl) {
      Disable (gbp->is_focus);
    }
    AlignObjects (ALIGN_CENTER, (HANDLE) g, (HANDLE) f, NULL);

    Hide (gbp->orgGrp [1]);

    gbp->orgGrp [2] = HiddenGroup (k, -1, 0, NULL);
    SetGroupSpacing (gbp->orgGrp [2], 10, 10);

    f = HiddenGroup (gbp->orgGrp [2], 2, 0, NULL);
    SetGroupSpacing (f, 10, 5);
    StaticPrompt (f, "Nuclear", 0, popupMenuHeight, programFont, 'l');
    gbp->gcode = PopupList (f, TRUE, NULL);
    PopulateGeneticCodePopup (gbp->gcode);
    SetValue (gbp->gcode, 1);
    StaticPrompt (f, "Mitochondrial", 0, popupMenuHeight, programFont, 'l');
    gbp->mgcode = PopupList (f, TRUE, NULL);
    PopulateGeneticCodePopup (gbp->mgcode);
    SetValue (gbp->mgcode, 1);

    Hide (gbp->orgGrp [2]);

    gbp->orgGrp [3] = HiddenGroup (k, -1, 0, NULL);
    SetGroupSpacing (gbp->orgGrp [3], 10, 10);

    f = HiddenGroup (gbp->orgGrp [3], 0, 2, NULL);
    y = StaticPrompt (f, "Taxonomic Lineage", 0, 0, programFont, 'c');
    gbp->lineage = ScrollText (f, 20, 3, programFont, TRUE, NULL);
    q = HiddenGroup (gbp->orgGrp [3], 2, 0, NULL);
    StaticPrompt (q, "Division", 0, dialogTextHeight, programFont, 'l');
    gbp->gbDiv = DialogText (q, "", 4, NULL);
    c = HiddenGroup (gbp->orgGrp [3], -1, 0, NULL);
    if (gfp != NULL && gfp->lookupTaxonomy != NULL) {
      ppt = StaticPrompt (c, "Looking up taxonomy will close this window.", 0, 0, programFont, 'c');
      b = PushButton (c, "Lookup Taxonomy", LookupTheTaxonomyProc);
      SetObjectExtra (b, gfp, NULL);
      AlignObjects (ALIGN_CENTER, (HANDLE) ppt, (HANDLE) b, NULL);
    }
    AlignObjects (ALIGN_CENTER, (HANDLE) f, (HANDLE) q, (HANDLE) c, NULL);

    Hide (gbp->orgGrp [3]);

    AlignObjects (ALIGN_CENTER, (HANDLE) tbs,
                  (HANDLE) gbp->orgGrp [0], (HANDLE) gbp->orgGrp [1],
                  (HANDLE) gbp->orgGrp [2], (HANDLE) gbp->orgGrp [3], NULL);

    pages [1] = HiddenGroup (m, -1, 0, NULL);
    SetGroupSpacing (pages [1], 10, 10);

    showUnstructMods = FALSE;
    hasSynonyms = FALSE;
    if (biop != NULL) {
      orp = biop->org;
      if (orp != NULL && orp->mod != NULL) {
        showUnstructMods = TRUE;
      }
      if (orp != NULL && orp->syn != NULL) {
        hasSynonyms = TRUE;
      }
    }

    if (showUnstructMods) {
      tbs = CreateFolderTabs (pages [1], modTabsUns, 0, 0, 0,
                              PROGRAM_FOLDER_TAB,
                              ChangeModSubPage, (Pointer) gbp);
    } else {
      tbs = CreateFolderTabs (pages [1], modTabs, 0, 0, 0,
                              PROGRAM_FOLDER_TAB,
                              ChangeModSubPage, (Pointer) gbp);
    }
    k = HiddenGroup (pages [1], 0, 0, NULL);

    gbp->modGrp [0] = HiddenGroup (k, -1, 0, NULL);
    SetGroupSpacing (gbp->modGrp [0], 10, 10);

    g = HiddenGroup (gbp->modGrp [0], -1, 0, NULL);
    SetGroupSpacing (g, 3, 10);
    gbp->subsrc_val_dlg = CreateSubSourceDialog (g, gbp->subsource_alists[0]);

    gbp->modGrp [1] = HiddenGroup (k, -1, 0, NULL);
    SetGroupSpacing (gbp->modGrp [1], 10, 10);

    g = HiddenGroup (gbp->modGrp [1], -1, 0, NULL);
    SetGroupSpacing (g, 3, 10);

    gbp->orgmod_val_dlg = CreateOrgModDialog (g, gbp->orgmod_alists[0], gbp->taxName);

    Hide (gbp->modGrp [1]);

    gbp->modGrp [2] = HiddenGroup (k, -1, 0, NULL);
    SetGroupSpacing (gbp->modGrp [2], 10, 10);

    g = HiddenGroup (gbp->modGrp [2], 2, 0, NULL);
    SetGroupSpacing (g, 3, 10);

    StaticPrompt (g, "Assigned Acronym", 0, stdLineHeight, programFont, 'l');
    gbp->gbacr = StaticPrompt (g, "", 15 * stdCharWidth, stdLineHeight, systemFont, 'l');
    StaticPrompt (g, "Assigned Anamorph", 0, stdLineHeight, programFont, 'l');
    gbp->gbana = StaticPrompt (g, "", 15 * stdCharWidth, stdLineHeight, systemFont, 'l');
    StaticPrompt (g, "Assigned Synonym", 0, stdLineHeight, programFont, 'l');
    gbp->gbsyn = StaticPrompt (g, "", 15 * stdCharWidth, stdLineHeight, systemFont, 'l');

    Hide (gbp->modGrp [2]);

    gbp->modGrp [3] = HiddenGroup (k, -1, 0, NULL);
    SetGroupSpacing (gbp->modGrp [3], 10, 10);

    if (showUnstructMods) {
      f3 = HiddenGroup (gbp->modGrp [3], 0, 2, NULL);
      StaticPrompt (f3, "Unstructured Modifiers", 0, 0, programFont, 'c');
      gbp->mod = CreateVisibleStringDialog (f3, 3, -1, 15);
    }

    Hide (gbp->modGrp [3]);

    AlignObjects (ALIGN_CENTER, (HANDLE) tbs,
                  (HANDLE) gbp->modGrp [0],
                  (HANDLE) gbp->modGrp [1],
                  (HANDLE) gbp->modGrp [2],
                  (HANDLE) gbp->modGrp [3], NULL);

    Hide (pages [1]);

    pages [2] = HiddenGroup (m, -1, 0, NULL);
    SetGroupSpacing (pages [2], 10, 10);

    tabs = miscTabs1;
    if (hasSynonyms) {
      tabs = miscTabs2;
    }
    tbs = CreateFolderTabs (pages [2], tabs, 0, 0, 0,
                            PROGRAM_FOLDER_TAB,
                            ChangeMiscSubPage, (Pointer) gbp);
    k = HiddenGroup (pages [2], 0, 0, NULL);

    for (z = 0; z < 2; z++) {
      gbp->miscGrp [z] = NULL;
    }
    z = 0;

    if (hasSynonyms) {
      gbp->miscGrp [z] = HiddenGroup (k, -1, 0, NULL);
      SetGroupSpacing (gbp->miscGrp [z], 10, 10);

      f1 = HiddenGroup (gbp->miscGrp [z], 0, 2, NULL);
      StaticPrompt (f1, "Synonyms", 0, 0, programFont, 'c');
      gbp->syn = CreateVisibleStringDialog (f1, 3, -1, 15);

      z++;
    }

    gbp->miscGrp [z] = HiddenGroup (k, -1, 0, NULL);
    SetGroupSpacing (gbp->miscGrp [z], 10, 10);

    f2 = HiddenGroup (gbp->miscGrp [z], -1, 0, NULL);
    SetGroupSpacing (f2, 10, 10);
    if (GetAppProperty ("ReadOnlyDbTags") == NULL) {
      just = 'c';
    } else {
      just = 'l';
      StaticPrompt (f2, "This page is read-only", 15 * stdCharWidth, 0, programFont, 'c');
    }
    t = HiddenGroup (f2, 2, 0, NULL);
    StaticPrompt (t, "Database", 7 * stdCharWidth, 0, programFont, just);
    StaticPrompt (t, "Object ID", 8 * stdCharWidth, 0, programFont, just);
    gbp->db = CreateDbtagDialog (f2, 3, -1, 7, 8);

    Hide (gbp->miscGrp [1]);

    AlignObjects (ALIGN_CENTER, (HANDLE) tbs,
                  (HANDLE) gbp->miscGrp [0],
                  (HANDLE) gbp->miscGrp [1], NULL);

    Hide (pages [2]);


    AlignObjects (ALIGN_CENTER, (HANDLE) pages [ORGANISM_PAGE],
                  (HANDLE) pages [MODIFIERS_PAGE],
                  (HANDLE) pages [MISCELLANEOUS_PAGE], NULL);
  }

  return (DialoG) p;
}

static void SetBioSourceImportExportItems (GenBioFormPtr gfp)

{
  IteM  exportItm;
  IteM  importItm;

  if (gfp != NULL) {
    importItm = FindFormMenuItem ((BaseFormPtr) gfp, VIB_MSG_IMPORT);
    exportItm = FindFormMenuItem ((BaseFormPtr) gfp, VIB_MSG_EXPORT);
    switch (gfp->currentPage) {
      case ORGANISM_PAGE :
        SafeSetTitle (importItm, "Import BioSource...");
        SafeSetTitle (exportItm, "Export BioSource...");
        SafeEnable (importItm);
        SafeEnable (exportItm);
        break;
      case MODIFIERS_PAGE :
        SafeSetTitle (importItm, "Import...");
        SafeSetTitle (exportItm, "Export...");
        SafeDisable (importItm);
        SafeDisable (exportItm);
        break;
      case MISCELLANEOUS_PAGE :
        SafeSetTitle (importItm, "Import...");
        SafeSetTitle (exportItm, "Export...");
        SafeDisable (importItm);
        SafeDisable (exportItm);
        break;
      case COMMON_PAGE :
        SafeSetTitle (importItm, "Import...");
        SafeSetTitle (exportItm, "Export...");
        SafeDisable (importItm);
        SafeDisable (exportItm);
        break;
      case LOCATION_PAGE :
        SafeSetTitle (importItm, "Import SeqLoc...");
        SafeSetTitle (exportItm, "Export SeqLoc...");
        SafeEnable (importItm);
        SafeEnable (exportItm);
        break;
      default :
        break;
    }
  }
}

static void ChangeBioSourcePage (VoidPtr data, Int2 newval, Int2 oldval)

{
  GenBioFormPtr  gfp;

  gfp = (GenBioFormPtr) data;
  if (gfp != NULL) {
    gfp->currentPage = newval;
    SafeHide (gfp->pages [oldval]);
    SafeShow (gfp->pages [newval]);
    switch (newval) {
      case ORGANISM_PAGE :
        break;
      case MODIFIERS_PAGE :
        break;
      case MISCELLANEOUS_PAGE :
        break;
      case COMMON_PAGE :
        break;
      case LOCATION_PAGE :
        SendMessageToDialog (gfp->location, VIB_MSG_ENTER);
        break;
      default :
        break;
    }
    SetBioSourceImportExportItems (gfp);
    Update ();
  }
}

static Boolean ImportBioSourceForm (ForM f, CharPtr filename)

{
  AsnIoPtr       aip;
  BioSourcePtr   biop;
  GenBioFormPtr  gfp;
  Char           path [PATH_MAX];

  path [0] = '\0';
  StringNCpy_0 (path, filename, sizeof (path));
  gfp = (GenBioFormPtr) GetObjectExtra (f);
  if (gfp != NULL) {
    switch (gfp->currentPage) {
      case ORGANISM_PAGE :
        if (path [0] != '\0' || GetInputFileName (path, sizeof (path), "", "TEXT")) {
          aip = AsnIoOpen (path, "r");
          if (aip != NULL) {
            biop = BioSourceAsnRead (aip, NULL);
            AsnIoClose (aip);
            if (biop != NULL) {
              PointerToDialog (gfp->data, (Pointer) biop);
              biop = BioSourceFree (biop);
              Update ();
              return TRUE;
            }
          }
        }
        break;
      case LOCATION_PAGE :
        return ImportDialog (gfp->location, filename);
      default :
        break;
    }
  }
  return FALSE;
}

static Boolean ExportBioSourceForm (ForM f, CharPtr filename)

{
  AsnIoPtr       aip;
  BioSourcePtr   biop;
  GenBioFormPtr  gfp;
  Char           path [PATH_MAX];
#ifdef WIN_MAC
  FILE           *fp;
#endif

  path [0] = '\0';
  StringNCpy_0 (path, filename, sizeof (path));
  gfp = (GenBioFormPtr) GetObjectExtra (f);
  if (gfp != NULL) {
    switch (gfp->currentPage) {
      case ORGANISM_PAGE :
        if (path [0] != '\0' || GetOutputFileName (path, sizeof (path), NULL)) {
#ifdef WIN_MAC
          fp = FileOpen (path, "r");
          if (fp != NULL) {
            FileClose (fp);
          } else {
            FileCreate (path, "TEXT", "ttxt");
          }
#endif
          aip = AsnIoOpen (path, "w");
          if (aip != NULL) {
            biop = DialogToPointer (gfp->data);
            BioSourceAsnWrite (biop, aip, NULL);
            AsnIoClose (aip);
            biop = BioSourceFree (biop);
            return TRUE;
          }
        }
        break;
      case LOCATION_PAGE :
        return ExportDialog (gfp->location, filename);
      default :
        break;
    }
  }
  return FALSE;
}

static CharPtr  biosourceDescFormTabs [] = {
  "Organism", "Modifiers", "Miscellaneous", NULL
};

static void BioSourceDescFormMessage (ForM f, Int2 mssg)

{
  GenBioFormPtr  gfp;

  gfp = (GenBioFormPtr) GetObjectExtra (f);
  if (gfp != NULL) {
    switch (mssg) {
      case VIB_MSG_INIT :
        SendMessageToDialog (gfp->data, VIB_MSG_INIT);
        break;
      case VIB_MSG_IMPORT :
        ImportBioSourceForm (f, NULL);
        break;
      case VIB_MSG_EXPORT :
        ExportBioSourceForm (f, NULL);
        break;
      case VIB_MSG_CLOSE :
        Remove (f);
        break;
      case VIB_MSG_CUT :
        StdCutTextProc (NULL);
        break;
      case VIB_MSG_COPY :
        StdCopyTextProc (NULL);
        break;
      case VIB_MSG_PASTE :
        StdPasteTextProc (NULL);
        break;
      case VIB_MSG_DELETE :
        StdDeleteTextProc (NULL);
        break;
      default :
        if (gfp->appmessage != NULL) {
          gfp->appmessage (f, mssg);
        }
        break;
    }
  }
}

static void BioSourceFormActivate (WindoW w)

{
  GenBioFormPtr  gfp;

  gfp = (GenBioFormPtr) GetObjectExtra (w);
  if (gfp != NULL) {
    if (gfp->activate != NULL) {
      gfp->activate (w);
    }
    SetBioSourceImportExportItems (gfp);
  }
}

static Boolean OkayToAcceptBioSource (ButtoN b)

{
  Boolean        abort;
  MsgAnswer      ans;
  OrgModPtr      curr;
  GenBioPagePtr  gbp;
  GenBioFormPtr  gfp;
  OrgModPtr      mod;
  CharPtr        str;
  ValNodePtr     err_list;

  gfp = (GenBioFormPtr) GetObjectExtra (b);
  if (gfp != NULL) {
    gbp = (GenBioPagePtr) GetObjectExtra (gfp->data);
    if (gbp != NULL) {
      err_list = TestDialog (gbp->orgmod_val_dlg);
      ValNodeLink (&err_list, TestDialog (gbp->subsrc_val_dlg));
      if (err_list != NULL)
      {
        if (ANS_CANCEL == Message (MSG_OKC, "You have selected values for modifiers, but not types - values without types will be discarded.  Continue?"))
        {
          err_list = ValNodeFreeData (err_list);
          return FALSE;
        }
        err_list = ValNodeFreeData (err_list);
      }
      mod = DialogToPointer (gbp->orgmod_val_dlg);
      if (mod == NULL) return TRUE;
      abort = TRUE;
      for (curr = mod; curr != NULL; curr = curr->next) {
        if (curr->subtype == 254) {
          abort = FALSE;
        }
      }
      OrgModFree (mod);
      if (abort) return TRUE;
      str = SaveStringFromText (gbp->taxName);
      ans = ANS_NO;
      if (StringICmp (str, gbp->origTaxName) != 0) {
        ans = Message (MSG_YNC, "Delete original name (necessary for correct lookup)?");
      }
      MemFree (str);
      if (ans == ANS_CANCEL) return FALSE;
      if (ans == ANS_YES) {
        gbp->stripOldName = TRUE;
      }
    }
  }
  return TRUE;
}

static void BioSourceDescFormAcceptButtonProc (ButtoN b)

{
  if (OkayToAcceptBioSource (b)) {
    StdAcceptFormButtonProc (b);
  }
}

extern ForM CreateBioSourceDescForm (Int2 left, Int2 top, Int2 width,
                                     Int2 height, CharPtr title, ValNodePtr sdp,
                                     SeqEntryPtr sep, FormActnFunc actproc,
                                     BioSourceEditProcsPtr bepp)

{
  ButtoN             b;
  BioSourcePtr       biop;
  GrouP              c;
  GrouP              g;
  GenBioFormPtr      gfp;
  StdEditorProcsPtr  sepp;
  WindoW             w;

  w = NULL;
  gfp = (GenBioFormPtr) MemNew (sizeof (GenBioForm));
  if (gfp != NULL) {
    w = FixedWindow (left, top, width, height, title, StdCloseWindowProc);
    SetObjectExtra (w, gfp, StdDescFormCleanupProc);
    gfp->form = (ForM) w;
    gfp->actproc = actproc;
    gfp->toform = NULL;
    gfp->fromform = NULL;
    gfp->testform = NULL;
    gfp->formmessage = BioSourceDescFormMessage;
    gfp->importform = ImportBioSourceForm;
    gfp->exportform = ExportBioSourceForm;

#ifndef WIN_MAC
    CreateStdEditorFormMenus (w);
#endif

    if (bepp == NULL) {
      bepp = (BioSourceEditProcsPtr) GetAppProperty ("BioSourcEditForm");
    }
    if (bepp != NULL) {
      gfp->lookupTaxonomy = bepp->lookupTaxonomy;
    }
    sepp = (StdEditorProcsPtr) GetAppProperty ("StdEditorForm");
    if (sepp != NULL) {
      gfp->activate = sepp->activateForm;
      gfp->appmessage = sepp->handleMessages;
    }
    SetActivate (w, BioSourceFormActivate);

    g = HiddenGroup (w, -1, 0, NULL);
    SetGroupSpacing (g, 3, 10);

    gfp->foldertabs = CreateFolderTabs (g, biosourceDescFormTabs, ORGANISM_PAGE,
                                        0, 0, SYSTEM_FOLDER_TAB,
                                        ChangeBioSourcePage, (Pointer) gfp);
    gfp->currentPage = ORGANISM_PAGE;

    biop = NULL;
    if (sdp != NULL && sdp->choice == Seq_descr_source) {
      biop = sdp->data.ptrvalue;
    }
    gfp->data = CreateBioSourceDialog (g, NULL, gfp->pages, biop, gfp, FALSE);

    AlignObjects (ALIGN_CENTER, (HANDLE) gfp->foldertabs, (HANDLE) gfp->data, NULL);

    c = HiddenGroup (w, 2, 0, NULL);
    b = PushButton (c, "Accept", BioSourceDescFormAcceptButtonProc);
    SetObjectExtra (b, gfp, NULL);
    PushButton (c, "Cancel", StdCancelButtonProc);
    AlignObjects (ALIGN_CENTER, (HANDLE) g, (HANDLE) c, NULL);
    RealizeWindow (w);

    SendMessageToDialog (gfp->data, VIB_MSG_INIT);
    Show (gfp->pages [gfp->currentPage]);
    SendMessageToDialog (gfp->data, VIB_MSG_ENTER);
    Update ();
  }
  return (ForM) w;
}

typedef struct biosourcedlg
{
  DIALOG_MESSAGE_BLOCK
  DialoG             foldertabs;
  GrouP              pages [NUM_PAGES];
  DialoG             location;
  Int2               currentPage;
  DialoG             gbp;
  
  LookupTaxonomyProc lookupTaxonomy; /* read in from the application properties */
} BioSourceDlgData, PNTR BioSourceDlgPtr;

static void ChangeBioSourceDialogPage (VoidPtr data, Int2 newval, Int2 oldval)

{
  BioSourceDlgPtr dlg;

  dlg = (BioSourceDlgPtr) data;
  if (dlg != NULL) {
    dlg->currentPage = newval;
    SafeHide (dlg->pages [oldval]);
    SafeShow (dlg->pages [newval]);
    switch (newval) {
      case ORGANISM_PAGE :
        break;
      case MODIFIERS_PAGE :
        break;
      case MISCELLANEOUS_PAGE :
        break;
      case COMMON_PAGE :
        break;
      case LOCATION_PAGE :
        SendMessageToDialog (dlg->location, VIB_MSG_ENTER);
        break;
      default :
        break;
    }
    Update ();
  }
}

static void BioSourceToDialog (DialoG d, Pointer userdata)
{
  BioSourceDlgPtr dlg;
  BioSourcePtr    biop;
  
  dlg = (BioSourceDlgPtr) GetObjectExtra (d);
  if (dlg == NULL)
  {
    return;
  }
  
  biop = (BioSourcePtr) userdata;
  PointerToDialog (dlg->gbp, biop);
}

static Pointer DialogToBioSource (DialoG d)
{
  BioSourceDlgPtr dlg;
  
  dlg = (BioSourceDlgPtr) GetObjectExtra (d);
  if (dlg == NULL)
  {
    return NULL;
  }
  
  return DialogToPointer (dlg->gbp);
}

static ValNodePtr TestBioSourceDialog (DialoG d)
{
  BioSourceDlgPtr dlg;
  
  dlg = (BioSourceDlgPtr) GetObjectExtra (d);
  if (dlg == NULL)
  {
    return NULL;
  }
  return TestDialog (dlg->gbp);
}

static void BioSourceDialogMessage (DialoG d, Int2 message)
{
  BioSourceDlgPtr dlg;
  
  dlg = (BioSourceDlgPtr) GetObjectExtra (d);
  if (dlg != NULL)
  {
    SendMessageToDialog (dlg->gbp, message);
  }    
}

static void DoTaxLookup (ButtoN b)
{
  BioSourceDlgPtr dlg;
  SeqEntryPtr     sep;
  BioseqSetPtr    bssp;
  SeqDescrPtr     sdp;
  BioSourcePtr    biop;
  Uint2           entityID;
 
  dlg = (BioSourceDlgPtr) GetObjectExtra (b);
  if (dlg == NULL)
  {
    return;
  }

  sep = SeqEntryNew ();
  if (sep == NULL)
  {
    return;
  }
  
  biop = (BioSourcePtr) DialogToPointer (dlg->dialog);
  bssp = BioseqSetNew ();
  
  if (biop == NULL || bssp == NULL)
  {
    biop = BioSourceFree (biop);
    bssp = BioseqSetFree (bssp);
    sep = SeqEntryFree (sep);
    return;
  }
  
  bssp->_class = BioseqseqSet_class_empty_set;
  
  sep->choice = 2;
  sep->data.ptrvalue = bssp;
  sdp = CreateNewDescriptor (sep, Seq_descr_source);
  sdp->data.ptrvalue = biop;
  entityID = ObjMgrGetEntityIDForChoice (sep);
  if (entityID > 0)
  {
    SeqMgrIndexFeatures (entityID, NULL);
    if ((dlg->lookupTaxonomy) (entityID))
    {
      /* find the new BioSource */
      sdp = bssp->descr;
      while (sdp != NULL && sdp->choice != Seq_descr_source)
      {
        sdp = sdp->next;
      }
      if (sdp != NULL)
      {
        /* put it in our dialog */
        PointerToDialog (dlg->dialog, sdp->data.ptrvalue);
      }
    }
  }                 
  SeqEntryFree (sep);
}

extern DialoG BioSourceDialog (GrouP parent)
{
  BioSourceDlgPtr       dlg;
  GrouP                 p;
  BioSourceEditProcsPtr bepp;
  ButtoN                b;
  
  dlg = (BioSourceDlgPtr) MemNew (sizeof (BioSourceDlgData));
  if (dlg == NULL)
  {
    return NULL;
  }
  
  p = HiddenGroup (parent, -1, 0, NULL);
  SetObjectExtra (p, dlg, StdCleanupExtraProc);
  SetGroupSpacing (p, 3, 10);
  dlg->dialog = (DialoG) p;
  dlg->todialog = BioSourceToDialog;
  dlg->fromdialog = DialogToBioSource;
  dlg->dialogmessage = BioSourceDialogMessage;
  dlg->testdialog = TestBioSourceDialog;

  dlg->foldertabs = CreateFolderTabs (p, biosourceDescFormTabs, ORGANISM_PAGE,
                                      0, 0, SYSTEM_FOLDER_TAB,
                                      ChangeBioSourceDialogPage, (Pointer) dlg);
  dlg->currentPage = ORGANISM_PAGE;

  dlg->gbp = CreateBioSourceDialog (p, NULL, dlg->pages, NULL, NULL, TRUE);
  
  bepp = (BioSourceEditProcsPtr) GetAppProperty ("BioSourcEditForm");
  
  if (bepp != NULL && bepp->lookupTaxonomy != NULL) {
    dlg->lookupTaxonomy = bepp->lookupTaxonomy;
    b = PushButton (p, "Look up Taxonomy", DoTaxLookup);
    SetObjectExtra (b, dlg, NULL);
  }
  else
  {
    b = NULL;
  }

  
  AlignObjects (ALIGN_CENTER, (HANDLE) dlg->foldertabs,
                              (HANDLE) dlg->gbp, 
                              (HANDLE) b, 
                              NULL);

  SendMessageToDialog (dlg->gbp, VIB_MSG_INIT);
  SendMessageToDialog (dlg->gbp, VIB_MSG_ENTER);

  return (DialoG) p;
}

static CharPtr  biosourceFeatFormTabs [] = {
  "Organism", "Modifiers", "Miscellaneous", "Properties", "Location", NULL
};

static void BioSourceFeatFormMessage (ForM f, Int2 mssg)

{
  GenBioFormPtr  gfp;

  gfp = (GenBioFormPtr) GetObjectExtra (f);
  if (gfp != NULL) {
    switch (mssg) {
      case VIB_MSG_INIT :
        SendMessageToDialog (gfp->data, VIB_MSG_INIT);
        StdInitFeatFormProc (f);
        break;
      case VIB_MSG_IMPORT :
        ImportBioSourceForm (f, NULL);
        break;
      case VIB_MSG_EXPORT :
        ExportBioSourceForm (f, NULL);
        break;
      case VIB_MSG_CLOSE :
        Remove (f);
        break;
      case VIB_MSG_CUT :
        StdCutTextProc (NULL);
        break;
      case VIB_MSG_COPY :
        StdCopyTextProc (NULL);
        break;
      case VIB_MSG_PASTE :
        StdPasteTextProc (NULL);
        break;
      case VIB_MSG_DELETE :
        StdDeleteTextProc (NULL);
        break;
      default :
        if (gfp->appmessage != NULL) {
          gfp->appmessage (f, mssg);
        }
        break;
    }
  }
}

static void BioSourceFeatFormAcceptButtonProc (ButtoN b)

{
  if (OkayToAcceptBioSource (b)) {
    StdFeatFormAcceptButtonProc (b);
  }
}

extern ForM CreateBioSourceFeatForm (Int2 left, Int2 top, Int2 width,
                                     Int2 height, CharPtr title, SeqFeatPtr sfp,
                                     SeqEntryPtr sep, FormActnFunc actproc,
                                     BioSourceEditProcsPtr bepp)

{
  ButtoN             b;
  BioSourcePtr       biop;
  GrouP              c;
  GrouP              g;
  GenBioFormPtr      gfp;
  GrouP              h;
  GrouP              s;
  StdEditorProcsPtr  sepp;
  WindoW             w;

  w = NULL;
  gfp = (GenBioFormPtr) MemNew (sizeof (GenBioForm));
  if (gfp != NULL) {
    w = FixedWindow (left, top, width, height, title, StdCloseWindowProc);
    SetObjectExtra (w, gfp, StdFeatFormCleanupProc);
    gfp->form = (ForM) w;
    gfp->actproc = actproc;
    gfp->toform = StdSeqFeatPtrToFeatFormProc;
    gfp->fromform = NULL;
    gfp->testform = NULL;
    gfp->formmessage = BioSourceFeatFormMessage;
    gfp->importform = ImportBioSourceForm;
    gfp->exportform = ExportBioSourceForm;

#ifndef WIN_MAC
    CreateStdEditorFormMenus (w);
#endif

    if (bepp == NULL) {
      bepp = (BioSourceEditProcsPtr) GetAppProperty ("BioSourcEditForm");
    }
    if (bepp != NULL) {
      gfp->lookupTaxonomy = bepp->lookupTaxonomy;
    }
    sepp = (StdEditorProcsPtr) GetAppProperty ("StdEditorForm");
    if (sepp != NULL) {
      gfp->activate = sepp->activateForm;
      gfp->appmessage = sepp->handleMessages;
    }
    SetActivate (w, BioSourceFormActivate);

    g = HiddenGroup (w, -1, 0, NULL);
    SetGroupSpacing (g, 3, 10);

    gfp->foldertabs = CreateFolderTabs (g, biosourceFeatFormTabs, ORGANISM_PAGE,
                                        0, 0, SYSTEM_FOLDER_TAB,
                                        ChangeBioSourcePage, (Pointer) gfp);
    gfp->currentPage = ORGANISM_PAGE;

    h = HiddenGroup (g, 0, 0, NULL);

    biop = NULL;
    if (sfp != NULL && sfp->data.choice == SEQFEAT_BIOSRC) {
      biop = sfp->data.value.ptrvalue;
    }
    gfp->data = CreateBioSourceDialog (h, NULL, gfp->pages, biop, gfp, TRUE);
 
    s = HiddenGroup (h, -1, 0, NULL);
    CreateCommonFeatureGroup (s, (FeatureFormPtr) gfp, sfp, FALSE, TRUE);
    gfp->pages [COMMON_PAGE] = s;
    Hide (gfp->pages [COMMON_PAGE]);

    s = HiddenGroup (h, -1, 0, NULL);
    gfp->location = CreateIntervalEditorDialogEx (s, NULL, 4, 2, sep, TRUE, TRUE,
                                                  TRUE, TRUE, FALSE,
                                                  (FeatureFormPtr) gfp,
                                                  StdFeatIntEdPartialCallback);
    gfp->pages [LOCATION_PAGE] = s;
    Hide (gfp->pages [LOCATION_PAGE]);

    AlignObjects (ALIGN_CENTER, (HANDLE) gfp->data,
                  (HANDLE) gfp->pages [COMMON_PAGE],
                  (HANDLE) gfp->pages [LOCATION_PAGE], NULL);
    AlignObjects (ALIGN_CENTER, (HANDLE) gfp->foldertabs, (HANDLE) h, NULL);

    c = HiddenGroup (w, 2, 0, NULL);
    b = PushButton (c, "Accept", BioSourceFeatFormAcceptButtonProc);
    SetObjectExtra (b, gfp, NULL);
    PushButton (c, "Cancel", StdCancelButtonProc);
    AlignObjects (ALIGN_CENTER, (HANDLE) g, (HANDLE) c, NULL);
    RealizeWindow (w);

    SendMessageToDialog (gfp->data, VIB_MSG_INIT);
    SendMessageToDialog (gfp->location, VIB_MSG_INIT);
    Show (gfp->pages [gfp->currentPage]);
    SendMessageToDialog (gfp->data, VIB_MSG_ENTER);
    Update ();
  }
  return (ForM) w;
}

extern Int2 LIBCALLBACK BioSourceGenFunc (Pointer data)

{
  GenBioFormPtr     gfp;
  HelpMessageFunc   helpfunc;
  Uint2             itemtype;
  OMProcControlPtr  ompcp;
  OMUserDataPtr     omudp;
  ObjMgrProcPtr     proc;
  ValNodePtr        sdp;
  SeqEntryPtr       sep;
  SeqFeatPtr        sfp;
  Uint2             subtype;
  WindoW            w;

  ompcp = (OMProcControlPtr) data;
  w = NULL;
  sdp = NULL;
  sfp = NULL;
  sep = NULL;
  itemtype = 0;
  subtype = 0;
  if (ompcp == NULL || ompcp->proc == NULL) return OM_MSG_RET_ERROR;
  proc = ompcp->proc;
  switch (ompcp->input_itemtype) {
    case OBJ_SEQFEAT :
      sfp = (SeqFeatPtr) ompcp->input_data;
      if (sfp != NULL && sfp->data.choice != SEQFEAT_BIOSRC) {
        return OM_MSG_RET_ERROR;
      }
      itemtype = OBJ_SEQFEAT;
      subtype = FEATDEF_BIOSRC;
      break;
    case OBJ_SEQDESC :
      sdp = (ValNodePtr) ompcp->input_data;
      if (sdp != NULL && sdp->choice != Seq_descr_source) {
        return OM_MSG_RET_ERROR;
      }
      itemtype = OBJ_SEQDESC;
      subtype = Seq_descr_source;
      break;
    case OBJ_BIOSEQ :
      break;
    case OBJ_BIOSEQSET :
      break;
    case 0 :
      break;
    default :
      return OM_MSG_RET_ERROR;
  }
  omudp = ItemAlreadyHasEditor (ompcp->input_entityID, ompcp->input_itemID,
                                ompcp->input_itemtype, ompcp->proc->procid);
  if (omudp != NULL) {
    gfp = (GenBioFormPtr) omudp->userdata.ptrvalue;
    if (gfp != NULL) {
      Select (gfp->form);
    }
    return OM_MSG_RET_DONE;
  }
  sep = GetTopSeqEntryForEntityID (ompcp->input_entityID);
  if (sfp != NULL) {
    w = (WindoW) CreateBioSourceFeatForm (-50, -33, -10, -10,
                                          "Organism Information", sfp, sep,
                                          StdFeatFormActnProc, NULL);
  } else if (sdp != NULL) {
    w = (WindoW) CreateBioSourceDescForm (-50, -33, -10, -10,
                                          "Organism Information", sdp, sep,
                                          StdDescFormActnProc, NULL);
  } else {
    itemtype = proc->inputtype;
    subtype = proc->subinputtype;
    if (itemtype == OBJ_SEQFEAT && subtype == FEATDEF_BIOSRC) {
      w = (WindoW) CreateBioSourceFeatForm (-50, -33, -10, -10,
                                            "Organism Information", sfp, sep,
                                            StdFeatFormActnProc, NULL);
    } else if (itemtype == OBJ_SEQDESC && subtype == Seq_descr_source) {
      w = (WindoW) CreateBioSourceDescForm (-50, -33, -10, -10,
                                            "Organism Information", sdp, sep,
                                            StdDescFormActnProc, NULL);
    } else {
      return OM_MSG_RET_ERROR;
    }
  }
  gfp = (GenBioFormPtr) GetObjectExtra (w);
  if (gfp != NULL) {
    gfp->input_entityID = ompcp->input_entityID;
    gfp->input_itemID = ompcp->input_itemID;
    gfp->input_itemtype = ompcp->input_itemtype;
    gfp->this_itemtype = itemtype;
    gfp->this_subtype = subtype;
    gfp->procid = ompcp->proc->procid;
    gfp->proctype = ompcp->proc->proctype;
    gfp->userkey = OMGetNextUserKey ();
    omudp = ObjMgrAddUserData (ompcp->input_entityID, ompcp->proc->procid,
                               OMPROC_EDIT, gfp->userkey);
    if (omudp != NULL) {
      omudp->userdata.ptrvalue = (Pointer) gfp;
      omudp->messagefunc = StdVibrantEditorMsgFunc;
    }
    SendMessageToForm (gfp->form, VIB_MSG_INIT);
    if (sdp != NULL) {
      PointerToDialog (gfp->data, (Pointer) sdp->data.ptrvalue);
      SetClosestParentIfDuplicating ((BaseFormPtr) gfp);
    } else if (sfp != NULL) {
      PointerToForm (gfp->form, (Pointer) sfp);
      SetClosestParentIfDuplicating ((BaseFormPtr) gfp);
    } else if (itemtype == OBJ_SEQFEAT) {
      SetNewFeatureDefaultInterval ((FeatureFormPtr) gfp);
    }
  }
  Show (w);
  Select (w);
  helpfunc = (HelpMessageFunc) GetAppProperty ("HelpMessageProc");
  if (helpfunc != NULL) {
    helpfunc ("Biological Source", NULL);
  }
  return OM_MSG_RET_DONE;
}


static Int4 ChooseNext (Nlm_QualNameAssocPtr PNTR choice_list, Int4 num_in_list)
{
  Int4 best_choice = -1, i;

  for (i = 0; i < num_in_list; i++) {
    if (choice_list[i] == NULL || choice_list[i]->name == NULL) {
      /* do nothing */
    } else if (best_choice == -1) {
      /* no previous choice */
      best_choice = i;
    } else if (StringCmp (choice_list[i]->name, choice_list[best_choice]->name) <= 0) {
      best_choice = i;
    }
  }
  return best_choice;
}

extern EnumFieldAssocPtr GetModifiersEnum (Boolean get_subsource, Boolean get_orgmod, Boolean get_discouraged, Boolean get_discontinued)
{
  Int4 best_choice, k;
  EnumFieldAssocPtr eap;
  CharPtr newname;
  Uint2   newval;
  Boolean inserted_note = FALSE;
  Int4    num_eap = 2; /* count terminators */
  Nlm_QualNameAssocPtr choice_list[6];
  Int4                 num_choices = 0;
  

  /* initialize choice_list */
  for (k = 0; k < 6; k++) {
    choice_list[k] = NULL;
  }

  /* need size of both because despite the fact that we need only one blank and one terminator, we are adding in two notes */
  if (get_orgmod) {
    for (k = 0; current_orgmod_subtype_alist[k].name != NULL; k++) { 
      num_eap++;
    }
    choice_list[num_choices ++] = current_orgmod_subtype_alist + 1;
  }
  if (get_subsource) {
    for (k = 0; current_subsource_subtype_alist[k].name != NULL; k++) { 
      num_eap++;
    }
    choice_list[num_choices ++] = current_subsource_subtype_alist + 1;
  }


  if (get_discouraged) {
    if (get_orgmod) {
      for (k = 0; discouraged_orgmod_subtype_alist[k].name != NULL; k++) {
        num_eap++;
      }
      choice_list[num_choices++] = discouraged_orgmod_subtype_alist;
    }

    if (get_subsource) {
      for (k = 0; discouraged_subsource_subtype_alist[k].name != NULL; k++) {
        num_eap++;
      }   
      choice_list[num_choices++] = discouraged_subsource_subtype_alist;
    }
  }

  if (get_discontinued) {
    if (get_orgmod) {
      for (k = 0; discontinued_orgmod_subtype_alist[k].name != NULL; k++) {
        num_eap++;
      }
      choice_list[num_choices++] = discontinued_orgmod_subtype_alist;
    }
    if (get_subsource) {
      for (k = 0; discontinued_subsource_subtype_alist[k].name != NULL; k++) {
        num_eap++;
      }
      choice_list[num_choices++] = discontinued_subsource_subtype_alist;
    }
  }

  eap = (EnumFieldAssocPtr) MemNew (sizeof (EnumFieldAssoc) * num_eap);

  eap[0].name = " ";
  eap[0].value = 0;
  k = 1;
  while ((best_choice = ChooseNext(choice_list, num_choices)) != -1) {
    newname = choice_list[best_choice]->name;
    newval = choice_list[best_choice]->value;
    if (get_subsource && get_orgmod && (best_choice == 1 || best_choice == 3 || best_choice == 5)) {
      /* add 1000 to let calling function distinguish between subsource and orgmod */
      newval += 1000;
    }
    choice_list[best_choice]++;
    /* only add in notes if we are getting both subsource and orgmod */
    if (!inserted_note && get_subsource && get_orgmod && StringCmp (newname, "Note") > 0) {
      eap[k].name = "Note -- OrgMod";
      eap[k].value = ORGMOD_other;
      k++;
      eap[k].name = "Note -- SubSource";
      eap[k].value = SUBSRC_other + 1000;
      k++;
      inserted_note = TRUE;
    }
    eap[k].name = newname;
    eap[k].value = newval;
    k++;
  }
  eap[k].name = NULL;
  eap[k].value = 0;
  return eap;
}

extern EnumFieldAssocPtr GetSubSourceAndOrgModEnum (Boolean get_discouraged, Boolean get_discontinued)
{
  return GetModifiersEnum (TRUE, TRUE, get_discouraged, get_discontinued);
}


