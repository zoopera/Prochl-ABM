!
! XXXXXXXXXXXXXXXX
! XXXXXXXXXXXXXXXX
! XXX IAM CODE XXX
! XXXXXXXXXXXXXXXX
! XXXXXXXXXXXXXXXX
!
! Program: IAM
! Version: 12.0.GCS.a
! Date: 6/3/2017
! Author: Ferdi L. Hellweger
!         Department of Civil and Environmental Engineering
!         Northeastern University
!         Boston
!         ferdi@coe.neu.edu
!
!
! -------------
! ------------- 
! --- notes ---
! -------------
! -------------
!
! flh, 29.08.2020
! changes for Prochlorococcus dR/dC evolution project
! 1. oLM = 6-8
!    note: for oLM 8 added a_cfp, but not implemented for (1) hotstart read/write or (2) updating genome
! 2. fixed error, reading from opposite strand
!
! -----------------
! --- externals ---
! -----------------
!
! Ziggurat: r4_uni, r4_nor, r4_nor_setup, shr3
! https://people.sc.fsu.edu/~jburkardt/c_src/ziggurat_openmp/ziggurat_openmp.html
!
! Netlib random: random_gamma
! http://www.netlib.org/random/random.f90
! note: This is used for gamma distribution. We are not using it in final
! simulations. A problem was encountered compiling with this function on
! the NEU Discovery cluster. Those functions are therefore commented out.
!
!
! --------------------
! --- genetic code ---
! --------------------
!
! Source: https://en.wikipedia.org/wiki/DNA_codon_table
!
! Ala/A	GCT, GCC, GCA, GCG
! Arg/R	CGT, CGC, CGA, CGG, AGA, AGG
! Asn/N	AAT, AAC
! Asp/D	GAT, GAC
! Cys/C	TGT, TGC
! Gln/Q	CAA, CAG
! Glu/E	GAA, GAG
! Gly/G	GGT, GGC, GGA, GGG
! His/H	CAT, CAC
! Ile/I	ATT, ATC, ATA
! Leu/L	TTA, TTG, CTT, CTC, CTA, CTG
! Lys/K	AAA, AAG
! Met/M	ATG
! Phe/F	TTT, TTC
! Pro/P	CCT, CCC, CCA, CCG
! Ser/S	TCT, TCC, TCA, TCG, AGT, AGC
! Thr/T	ACT, ACC, ACA, ACG
! Trp/W	TGG
! Tyr/Y	TAT, TAC
! Val/V	GTT, GTC, GTA, GTG
! START	ATG	
! STOP	TAA, TGA, TAG
!
! integer equivalents, MW, C, N
!
! Source: https://en.wikipedia.org/wiki/Amino_acid
!
! Ala/A	1, 89, 3, 1
! Arg/R	2, 174, 6, 4
! Asn/N	3, 132, 4, 2
! Asp/D	4, 133, 4, 1
! Cys/C	5, 121, 3, 1
! Gln/Q	6, 147, 5, 1
! Glu/E	7, 146, 5, 2
! Gly/G	8, 75, 2, 1
! His/H	9, 155, 6, 3
! Ile/I	10, 131, 6, 1
! Leu/L	11, 131, 6, 1
! Lys/K	12, 146, 6, 2
! Met/M	13, 149, 5, 1
! Phe/F	14, 165, 9, 1
! Pro/P	15, 115, 5, 1
! Ser/S	16, 105, 3, 1
! Thr/T	17, 119, 4, 1
! Trp/W	18, 204, 11, 2
! Tyr/Y	19, 181, 9, 1
! Val/V	20, 117, 5, 1
! START	21
! STOP	22
!
! A  1
! T  2
! C  3
! G  4
!
!
! ----------------
! --- hotstart ---
! ----------------
!
! The present version supports hotstarting (checkpointing). This requires that
! the base genome does not change. i.e. all changes are stored at the agent
! level and no changes are flushed to the genome library. The model reads the
! general input (IAM_I_g.txt), not the hotstart genome (IAM_H_DNA.txt).
!

program IAM

!
! ---------------
! ---------------
! --- options ---
! ---------------
! ---------------
!
use omp_lib
!use random

implicit none

!
! -------------------------
! -------------------------
! --- declare variables ---
! -------------------------
! ------------------------- 


integer(4) cmt, crt

!
! ------------------
! --- dimensions ---
! ------------------
!
integer(4) dimna ! max. no. agents
integer(4) dimna_sp ! max. no. agents in subpopulation
integer(4) dimnc ! max. no of genome changes before update
integer(4) dimng ! max. no of genomes
integer(4) dimnnt ! max. no. bp in DNA
integer(4) dimnp ! max. no of proteins
integer(4) dimnsp ! no. subpopulations
integer(4) dimntb ! max. no. time variable inputs
integer(4) dimns ! no. values in sarkar lookup
integer(4) dimnk ! no. values in kimura lookup

!
! -------------------
! --- run control ---
! -------------------
!
! --- time ---
!
real(8) t ! time
real(8), allocatable :: t_sp(:)
real(8) dt ! time step
real(8) tend ! end time

!
! --- cluster output ---
!
integer(4) cTime ! cluster time
integer(4) cComp ! cluster compute node
integer(4) cIter ! cluster iteration
character(18) cString

!
! --- hotstart ---
!
integer(4) oHOT ! 0 = none, 1 = write, 2 = read, 3 = both
integer(4) oHOTs ! adjust sr from hotstart pop? 0 = no, 1 = yes
real(8) th ! total time
integer(4) iHOT ! hotstart counter
real(8) tHOT
integer(4) natHOT
real(8) SRTHOT
integer(4) npCHK ! number of checks
integer(4) icomp ! statistics for hotstart population by comp
integer(4) ncomp
integer(4) nat_comp(100)
real(8) SRT_comp(100)
integer(4), allocatable :: nat_comp_sp(:,:)
real(8), allocatable :: SRT_comp_sp(:,:)
!integer(4) ncompfrom
!integer(4) ncompto
!integer(4), allocatable :: icomptofrom(:,:)

!
! --- error handling --- 
!
integer(4), allocatable :: err_sp(:) ! error code
integer(4) ierr ! error flag

!
! --- HANS ---
!
integer(4) oHANS ! use HANS?, 0 = no, 1 = yes
integer(4) oHANSu ! distribution of mutations, 0 = diverse, 1 = uniform
integer(4) oHANSd ! use approaximation for dilution?, 0 = new draw, 1 = proportional
integer(4) oHANSs ! calculation method for number of mutants at end of time step:
! 1 = Luria & Delbrueck (1943)
! 2 = Sarkar et al (1992)
! 3 = Sarkar et al (1992) with preset & lookup
integer(4) oHANSr ! round sr of mother cell?, 0 = no, 1 = yes
integer(4) sarkar ! sarkar function
integer(4) isarkar ! sarkar flag, 0 = need to calc preset
integer(4), allocatable :: sarkarm(:) ! sarkar preset, mutation
integer(4), allocatable :: sarkarr(:) ! sarkar preset, recombination
real(8) sarkt, sarkn ! statistics
real(8) sarksr0

!
! --- iteration ---
!
integer(4) oIter ! iterate?, 0 = no, 1 = yes
integer(4) nIter ! no. iterations
integer(4) iIter ! iteration index

!
! --- random number ---
!
real(8) RSi, RSj ! random number seed inputs
integer(4) idum, jdum ! random number seeds
integer(4), allocatable :: idum_sp(:)
real(4) r ! random number
real(4) r4_uni ! ziggurat variables ...
real(4) r4_nor
real(4) fn(128)
integer(4) kn(128)
real(4) wn(128)
logical rfirst ! gamma distribution
integer(4) pcheck ! probability check function

!
! --- code performance --- 
!
real(8) cp_t1, cp_delta
real(8) cp_init_t1, cp_init_delta
real(8) cp_iter_t1, cp_iter_delta
real(8) cp_fix_t1, cp_fix_delta
real(8) cp_genomes_t1, cp_genomes_delta
real(8) cp_output_t1, cp_output_delta
real(8), allocatable :: cp_ENS_t1_sp(:), cp_ENS_delta_sp(:)
real(8) cp_ENS_delta
real(8), allocatable :: cp_HANS1_t1_sp(:), cp_HANS1_delta_sp(:)
real(8) cp_HANS1_delta
real(8), allocatable :: cp_HANS2_t1_sp(:), cp_HANS2_delta_sp(:)
real(8) cp_HANS2_delta
!real(8), allocatable :: cp_HANS2a_t1_sp(:), cp_HANS2a_delta_sp(:)
!real(8) cp_HANS2a_delta
!real(8), allocatable :: cp_HANS2b_t1_sp(:), cp_HANS2b_delta_sp(:)
!real(8) cp_HANS2b_delta
!real(8), allocatable :: cp_HANS2c_t1_sp(:), cp_HANS2c_delta_sp(:)
!real(8) cp_HANS2c_delta
real(8) cp_mix1_t1, cp_mix1_delta
real(8) cp_mix2_t1, cp_mix2_delta
real(8) cp_partot_t1, cp_partot_delta
real(8) cp_parwait_t1, cp_parwait_delta
real(8), allocatable :: cp_mr1_t1_sp(:), cp_mr1_delta_sp(:)
real(8) cp_mr1_delta
real(8), allocatable :: cp_mr2_t1_sp(:), cp_mr2_delta_sp(:)
real(8) cp_mr2_delta
!real(8), allocatable :: cp_mr2a_t1_sp(:), cp_mr2a_delta_sp(:)
!real(8) cp_mr2a_delta
!real(8), allocatable :: cp_mr2b_t1_sp(:), cp_mr2b_delta_sp(:)
!real(8) cp_mr2b_delta
!real(8), allocatable :: cp_mr2c_t1_sp(:), cp_mr2c_delta_sp(:)
!real(8) cp_mr2c_delta
real(8) cp_other_t1, cp_other_delta

!
! -------------
! --- cells ---
! -------------
!
integer(4) oSlim ! omit nonessential state variables?, 0 = no, 1 = yes
! note: not using this option, not tested fully (e.g. see POP print)

integer(4) ia1, ia2, ia3, ia4 ! agent indices

integer(4) isp1, isp2, isp3, isp4 ! subpopulation indices
integer(4) ntd ! number of threads
integer(4) nspe ! number of empty sp

integer(4) nat, naa, naf ! total, active, free agents
integer(4) nat_sp_mx
integer(4), allocatable :: nat_sp(:,:)
integer(4), allocatable :: naa_sp(:)
integer(4), allocatable :: naf_sp(:)

integer(1), allocatable :: a_sa(:,:,:) ! active
real*4, allocatable :: a_sr(:,:,:) ! super-individual up-scaling number
real*4, allocatable :: a_srm(:,:,:) ! temporary mutation counter/flag
real*4, allocatable :: a_srr(:,:,:) ! temporary recombination counter/flag
integer(4), allocatable :: a_iar(:,:,:) ! NOT USED
integer(4), allocatable :: a_gid(:,:,:) ! genome id 
integer(4), allocatable :: a_cn(:,:,:) ! change: number
integer(4), allocatable :: a_ck(:,:,:,:) ! change: position
integer(1), allocatable :: a_cnt(:,:,:,:) ! change: nt
real(4), allocatable :: a_cfp(:,:,:,:) ! change: fp
real*4, allocatable :: a_ntnk(:,:,:,:) ! no. of nucleotides in DNA by ATCG
integer(4), allocatable :: a_naa(:,:,:,:) ! no. amino acids
integer(4), allocatable :: a_naax(:,:,:,:,:) ! no. amino acid changes
real(8), allocatable :: a_uEg(:,:,:) ! mutation rate
real(8), allocatable :: a_q0DNAC(:,:,:)
real(8), allocatable :: a_q0DNAN(:,:,:)
real(8), allocatable :: a_q0DNAP(:,:,:)
real*4, allocatable :: a_q0aaC(:,:,:)
real*4, allocatable :: a_q0aaN(:,:,:)
real*4, allocatable :: a_VmaxC(:,:,:) ! max. C uptake rate
real*4, allocatable :: a_VmaxN(:,:,:) ! max. N uptake rate
real*4, allocatable :: a_VmaxP(:,:,:) ! max. P uptake rate

integer(4), allocatable :: a_sn(:,:,:) ! serial no.
integer(4), allocatable :: a_rsn(:,:,:) ! root serial no.
integer(4), allocatable :: a_nd(:,:,:) ! no. divisions
real*4, allocatable :: a_ntnSk(:,:,:,:) ! no. of S nucleotides in DNA by ATCG
real*4, allocatable :: a_ntnNk(:,:,:,:) ! no. of N nucleotides in DNA by ATCG
real*4, allocatable :: a_ntnXk(:,:,:,:) ! no. of X nucleotides in DNA by ATCG
integer(4), allocatable :: a_nmt(:,:,:) ! no. mutations, total
integer(4), allocatable :: a_nmS(:,:,:) ! no. mutations, synonymous
integer(4), allocatable :: a_nmN(:,:,:) ! no. mutations, nonsynonymous
integer(4), allocatable :: a_nrt(:,:,:) ! no. recombinations, total
integer(4), allocatable :: a_nrS(:,:,:) ! no. recombinations, synonymous
integer(4), allocatable :: a_nrN(:,:,:) ! no. recombinations, nonsynonymous
integer(1), allocatable :: a_cXSN(:,:,:,:) ! change: none (0), S (1) or N (2)
integer(1), allocatable :: a_cmr(:,:,:,:) ! change: none/undone (0), m (1) or r (2)
integer(4), allocatable :: a_nm(:,:,:,:,:) ! no. mutations, by nt
integer(4), allocatable :: a_nr(:,:,:,:,:) ! no. recombinations, by nt
real(8), allocatable :: a_nntt(:,:,:,:) ! no. nucleotides in mRNA
real(8), allocatable :: a_VC(:,:,:) ! C uptake rate
real(8), allocatable :: a_VN(:,:,:) ! N uptake rate
real(8), allocatable :: a_VP(:,:,:) ! P uptake rate
real(8), allocatable :: a_qC(:,:,:) ! C quota
real(8), allocatable :: a_qN(:,:,:) ! N quota
real(8), allocatable :: a_qP(:,:,:) ! P quota
real(8), allocatable :: a_CLimit(:,:,:) ! degree of C limitation (qC/q0C)
real(8), allocatable :: a_NLimit(:,:,:) ! degree of N limitation (qN/q0N)
real(8), allocatable :: a_PLimit(:,:,:) ! degree of P limitation (qP/q0P)
real(8), allocatable :: a_tb(:,:,:) ! birth time
real(8), allocatable :: a_tg(:,:,:) ! generation time

integer(4), allocatable :: iaf_sp(:,:) ! free agent index

integer(4), allocatable :: nar_sp(:) ! number of recombinations
integer(4), allocatable :: iad_sp(:,:) ! recombination donor index

!
! --- population size ---
!
real(8) SRT ! SR total
real(8), allocatable :: SRT_sp(:)
real(8) SRCDF ! SR CDF
real(8) SRCDFd ! SR CDF step
real(8) SRCDFs ! SR CDF start
real(8) SRCDFr ! SR CDF current
real(8) SRmin(100)
real(8) SRmax(100)
integer(4) iamax(100)
integer(4) iamin(100)
integer(4) imx1, imx2, nmx
real(8) fkill ! cells killed at mutation and recombination

!
! --- serial number ---
!
integer(4) gsn ! global serial number
! note: not coded to provide unique sn for multiple sp
integer(4), allocatable :: gsn_sp(:)
integer(4) gsnmax ! global serial number max

!
! --- initial population ---
!
integer(4) nic
integer(4) nic_sp
real(8) Popic
real(8), allocatable :: Popic_sp(:)

!
! -------------------------
! --- genome & proteins ---
! -------------------------
!
integer(4) oStartDNA ! 1 = diverse, random, 2 = uniform, random, 3 = uniform, all A, 4 = same, read

integer(4) ig1, ig2 ! genome indices
integer(4) ngt, nga ! total, active genomes

integer(4) ic1, ic2, ic3, icn1, icn2 ! change indices
integer(4) ncg_sp_mx ! number of changes
integer(4) ncg_sp_prev
integer(4), allocatable :: ncg_sp(:)
integer(4) icg, jcg ! genome updating flags

integer(4) iXSN ! change: -1 or -2 = none, 0 = X, 1 = S, 2 = N

integer(4) nt1, nt2, nt11, nt21, nt12, nt22 ! nucleotide
integer(4) nta(3), ntb(3), nta_o(3) ! 3 nucleotides
integer(4) nt5(5) ! 5 nucleotides

integer(4) int1, int2, int3, int4, int5 ! nucleotide indices
integer(4) nnt ! total nucleotides

real(8) bS(6), bN(6), bX(6), bSX(6), bT(6)
real(8) bSxu(6), bNxu(6), bXxu(6), bSXxu(6), bTxu(6)

integer(4) cod1, cod2, cod1_o ! codon

integer(4) gc(444) ! genetic code

integer(4) iaa1, iaa2, iaa1_o ! amino acid indices

integer(4) icod1 ! codon index

integer(4) ip1, ip2 ! protein indices
integer(4) ipw ! watch protein index

integer(4) pntstartx ! protein nucleotide stepping
integer(4) pntstopx
integer(4) pntstepx

integer(1), allocatable :: g_nt(:,:) ! nt
integer(1), allocatable :: g_nt_o(:,:) ! nt orig
integer(1), allocatable :: g_ntfS(:,:) ! nt fold synonymous
real(8), allocatable :: g_ntn(:) ! no. nt total
real(8), allocatable :: g_ntnS(:) ! no. nt synonymous
real(8), allocatable :: g_ntnN(:) ! no. nt nonsynonymous
real(8), allocatable :: g_ntnX(:) ! no. nt noncoding
real(8), allocatable :: g_ntnk(:,:) ! no. nt total by ATCG
real(8), allocatable :: g_ntnSk(:,:) ! no. nt synonymous by ATCG
real(8), allocatable :: g_ntnNk(:,:) ! no. nt nonsynonymous by ATCG
real(8), allocatable :: g_ntnXk(:,:) ! no. nt noncoding by ATCG
real(8), allocatable :: g_ntnij(:,:,:) ! no. nt total by ATCG by ATCG
real(8), allocatable :: g_ntnSij(:,:,:) ! no. nt synonymous by ATCG by ATCG
real(8), allocatable :: g_ntnNij(:,:,:) ! no. nt nonsynonymous by ATCG by ATCG
real(8), allocatable :: g_ntnXij(:,:,:) ! no. nt noncoding by ATCG by ATCG
real(8), allocatable :: g_nntt(:,:) ! no. nucleotides in mRNA
real(8), allocatable :: g_nnttt(:) ! total no. nucleotides in mRNA
real(8), allocatable :: g_nnttw(:,:) ! no. nucleotides in mRNA
real(8), allocatable :: g_nntttw(:) ! total no. nucleotides in mRNA
integer(4), allocatable :: g_naa(:,:) ! no. amino acids
real(8), allocatable :: g_uEg(:) 
integer(4), allocatable :: g_ip(:,:) ! protein id for this nt
integer(1), allocatable :: g_ntu1(:,:) ! nt updated for output and recombination
integer(1), allocatable :: g_ntu2(:,:)
integer(1), allocatable :: g_cmru1(:,:)
integer(1), allocatable :: g_cmru2(:,:)
integer(4), allocatable :: g_cic1(:,:) ! change no.
!integer(4), allocatable :: g_cic2(:,:)
real(4), allocatable :: g_cfp1(:,:) ! fp
real(4), allocatable :: g_cfp2(:,:)
integer(4), allocatable :: g_ncf(:) ! current fixed at cell level
integer(4), allocatable :: g_ncSf(:)
integer(4), allocatable :: g_ncNf(:)
integer(4), allocatable :: g_ncmf(:)
integer(4), allocatable :: g_ncrf(:)
integer(4) g_ncft
integer(4), allocatable :: g_ncff(:) ! flushed to genome
integer(4), allocatable :: g_ncSff(:)
integer(4), allocatable :: g_ncNff(:)
integer(4), allocatable :: g_ncmff(:)
integer(4), allocatable :: g_ncrff(:)
integer(1), allocatable :: g_cmr(:,:)
real(8), allocatable :: g_q0DNAC(:)
real(8), allocatable :: g_q0DNAN(:)
real(8), allocatable :: g_q0DNAP(:)
real(8), allocatable :: g_q0aaC(:)
real(8), allocatable :: g_q0aaN(:)
real(8), allocatable :: g_Vmax0C(:)
real(8), allocatable :: g_Vmax0N(:)
real(8), allocatable :: g_Vmax0P(:)
integer(4), allocatable :: g_pntstart(:,:) ! protein start
integer(4), allocatable :: g_pntstop(:,:) ! protein stop
integer(1), allocatable :: g_pother(:,:) ! protein encoded on other chain
real(8), allocatable :: g_pgamma(:,:) ! protein relative expression
real(8), allocatable :: g_pntn(:,:) ! protein no. nt total
real(8) pntnt
real(8), allocatable :: g_pntnS(:,:) ! protein no. nt synonymous
real(8) pntnSt
real(8), allocatable :: g_pntnN(:,:) ! protein no. nt nonsynonymous
real(8) pntnNt
real(8), allocatable :: g_pnmS(:,:) ! protein no. mutations, synonymous
real(8), allocatable :: g_pnmN(:,:) ! protein no. mutations, nonsynonymous
real(8), allocatable :: g_pnrS(:,:) ! protein no. recombinations, synonymous
real(8), allocatable :: g_pnrN(:,:) ! protein no. recombinations, nonsynonymous
integer(4), allocatable :: g_np(:) ! number of proteins
integer(4), allocatable :: g_naat(:) ! total no. aa coded

integer(4), allocatable :: ifx(:) ! fixed changes
integer(4), allocatable :: fx_cn(:)
integer(4), allocatable :: fx_ck(:,:)
integer(1), allocatable :: fx_cnt(:,:)
integer(1), allocatable :: fx_cXSN(:,:)
integer(1), allocatable :: fx_cmr(:,:)

integer(4), allocatable :: ifx_sp(:,:)
integer(4), allocatable :: fx_cn_sp(:,:)
integer(4), allocatable :: fx_ck_sp(:,:,:)
integer(1), allocatable :: fx_cnt_sp(:,:,:)
integer(1), allocatable :: fx_cXSN_sp(:,:,:)
integer(1), allocatable :: fx_cmr_sp(:,:,:)

integer(1) int_opo(1:4) ! nucleotide opposite strand

!
! --- reading & writing ----
!
integer(4) iger ! genetic element read index
integer(4) nger ! number of genetic elements to read
character(100) shr ! genome header read string
character(1), allocatable :: sntr(:) ! basepair read string
! note: multiple genetic elements read are stored as a single geneome

character(1) int_snt(1:4) ! nucleotide integer - string equivalents

integer(4) pntstartr ! for reading
integer(4) pntstopr
integer(4) potherr
real(8) pgammar

character(1) int_saa(1:22)
character(3) int_scod(1:64)
integer(4) int_icod(1:1000)

!
! --- synthetic ----
!
real(8) ficA ! initial fraction A
real(8) ficT ! initial fraction T
real(8) ficC ! initial fraction C
real(8) ficG ! initial fraction G

integer(4) pspacing ! protein spacing
integer(4) pminspacing
integer(4) pmaxspacing
integer(4) plength ! protein length
integer(4) pminlength
integer(4) pmaxlength

!
! --- perturbation analysis ----
!
real(8) pfAT ! initial fraction AT change to GC
real(8) pfGC ! initial fraction GC change to AT

!
! --- statistics ---
!
real(8) potherf ! protein other fraction
real(8) pgammat ! protein expression
integer(4) naaz(20) ! aa counter
integer(4) naazw(20) ! aa counter for watch protein
integer(4) naaxw
integer(4) ncodz(64) ! cod counter

!
! --------------------------
! --- reactor simulation ---
! --------------------------
!
! --- extracellular concentrations ---
!
integer(4) oFixExt ! fix extracellular concentration, 0 = no, 1 = constant, 2 = time variable
integer(4) oImpExt ! implicit integration for extracellular, 0 = no, 1 = yes

real(8) SC0 ! initial concentration
real(8) SN0
real(8) SP0

real(8) SinC, SinCB ! inflow concentration
real(8) SinN, SinNB
real(8) SinP, SinPB

real(8), dimension(:), allocatable :: SC_sp ! current concentration
real(8), dimension(:), allocatable :: SN_sp
real(8), dimension(:), allocatable :: SP_sp

real(8), dimension(:), allocatable :: RSC_sp ! derivative
real(8), dimension(:), allocatable :: RSN_sp
real(8), dimension(:), allocatable :: RSP_sp

!
! --- volume ---
!
real(8) VOL ! volume
real(8), dimension(:), allocatable :: VOL_sp

!
! --- flow rate ---
!
real(8) Q, QB ! inflow rate
real(8), dimension(:), allocatable :: Q_sp

!
! --- dilution ---
!
integer(4) opD ! dilution output option (0 = none, 1 = after D, 2 = before D)
real(8) dtD ! dilution time step
real(8) fD ! dilution fraction, fraction kept
integer(4) ipD ! dilution flag (1 = dilution just happened)
real(8) tDNext ! next dilution time
real(8), allocatable :: tDNext_sp(:)
real(8) bucci ! bucci poisson function
real(4) PD ! probability of dilution
real(8) PS ! probability of survival
real(8) nx ! number before dilution
real(8) hx ! expected number survival
real(8) rx ! actual number survival
integer(4) ix ! actual number survival
real(8) fx ! survival fraction
real(8) srx, srmx, srrx
integer(4) sax

!
! --- mixing ---
!
integer(4) oMix ! 1 = by SR, 2 = by agent
real(8) dtmix ! subpulation mix time step
integer(4) nspmix ! number of destination subpopulations to mix to
integer(4) ntx ! number of steps between mixing
integer(4) itx ! step index

!
! ------------------------------------------------------
! --- nutrient uptake and intracellular mass balance ---
! ------------------------------------------------------
!
real(8) Vmax0C
real(8) Vmax0N
real(8) Vmax0P

real(8), allocatable :: VmaxCa_sp(:) ! average maximum uptake rates
real(8), allocatable :: VmaxNa_sp(:)
real(8), allocatable :: VmaxPa_sp(:)

real(8) SUMVmax ! sum of Vmax

!
! --- effect of extracellular nutrient ---
!
integer(4) oLN ! limit uptake by extracellular nutrient, 0 = no, 1 = yes
real(8) KmC
real(8) KmN
real(8) KmP
real(8) LNC ! limitation factor, nutrient C
real(8) LNN ! limitation factor, nutrient N
real(8) LNP ! limitation factor, nutrient P

!
! --- effect of population ---
!
integer(4) oLP ! limit uptake by population, 0 = no, 1 = yes
real(8) K ! carrying capacity
real(8) P ! population size
real(8) LP ! limitation factor, population

!
! --- effect of mutations ---
!
integer(4) oLM
! limit uptake based on mutations?
! 0 = no
! 1 = constant
! 2 = exponential
! 3 = gamma
! 4 = constant & fraction neutral
! 5 = aa dissimilarity
integer(4) oKimura ! 1 = direct, 2 = lookup
real(8) fpn
real(8) fps ! fp for synonymous sites
real(8) fpx ! fp for non-coding sites
real(8) fpi
real(4) fp1, fp2
real(8) rnn 
real(8) rns ! rn for synonymous sites
real(8) rnx ! fp for non-coding sites
real(8) rni
real(8) si
integer(4) oFTSS
! mutations to/from START/STOP?
! 0 = OK
! 1 = kill
real(4) fpbeta ! kimura
real(8) spave
integer(4), allocatable :: kimura(:)

real(8) Mdc
real(8) Mfc
real(8) MdXX(22,22)

real(8) fdel

!
! --- other effect of GC ---
!
integer(4) oLoth ! 0 = no, 1 = yes
integer(4) oLothIC ! 1 = use init. genome for fGCt0, 2 = read in fGCt0
real(8) sGCt
real(8) fGCt
real(8) fGCt0
real(8) Loth ! limitation factor, other

!
! --- correction at division ---
!
integer(4) oCAD
real(8) dtad

!
! --- growth rate ---
!
real(8) kg, kgC, kgN, kgP ! kg actual, based on C, N, P
real(8) sr0 ! SR before growth
real(8) xsr ! SR factor

!
! ----------------
! --- division ---
!-----------------
!
real(8) fdq0 ! factor on change of q0
integer(4) oFixCells ! fix no. cells, 0 = no, 1 = constant, 2 = after dilution
integer(4) oLQ ! adjust q0 based on mutations?
! 0 = no
! 11 = DNA: change, use; aa: no change, no use
! 12 = DNA: change, no use; aa: no change, no use
! 21 = DNA: no change, no use; aa: change, use 
! 22 = DNA: no change, no use; aa: change, no use
! 31 (or 1) = DNA: change, use; aa: change, use
! 32 = DNA: change, no use; aa: change, no use
! 33 = DNA: change, use; aa: change, no use
! 34 = DNA: change, no use; aa: change, use
integer(4) oKillDaughter ! kill daughter cells?, 0 = no, 1 = yes

real(8) kreft_p ! kreft function 
real(8) SF ! split fraction
real(8) mqC ! mother qC
real(8) mqN ! mother qN
real(8) mqP ! mother qC
real(8) RCVmm ! RCV mean for mass
real(8) RCVmcv ! RCV CV for mass
real(8) RCVmcvl ! RCV CV limit for mass
real(8) RCVmll ! RCV low limit for mass
real(8) RCVmhl ! RCV high limit for mass

!
! --- nutrient quotas ---
!
real(8) NA ! Avacado constant

real(8) rGCC ! no. of C atoms in GC bp
real(8) rATC ! no. of C atoms in AT bp
real(8) rGCN ! no. of N atoms in GC bp
real(8) rATN ! no. of N atoms in AT bp
real(8) rGCP ! no. of P atoms in GC bp
real(8) rATP ! no. of P atoms in AT bp

!integer(4) MWaa(20) ! not used
real(8) raaC(20)
real(8) raaN(20)

integer(4) aaC
integer(4) aaN

real(8) q0aaCw
real(8) q0aaNw

real(8) n0DNA

real(8) m0aa

real(8) q0othC
real(8) q0othN
real(8) q0othP

real(8) qmaxC
real(8) qmaxN
real(8) qmaxP

real(8) fqmaxC
real(8) fqmaxN
real(8) fqmaxP

real(8) q0C
real(8) q0N
real(8) q0P

real(8) qaveC
real(8) qaveN
real(8) qaveP

!
! --------------------------------
! --- mutation & recombination ---
! --------------------------------
!
! --- mutation ---
!
real(8) fu ! factor on mutation rate (permanent)
real(8) uEg ! mutation error rates
real(8) uEgmin
real(8) uEgmax
real(8) ubij(4,4)
real(8) uEbij(4,4)
real(8) uEbk(4)
integer(4) cint(4,3) ! mutation destinations
real(4) P1 ! probability
integer(4) kmut ! mutation flag
real(8) fM ! factor on mutation rate (temporary)
integer(4) xM ! multiplier on mutations
integer(4) xMuse
integer(4) ixmr
real(8) fmKill ! fraction of mutants to kill (debug use)

!
! --- recombination ---
!
integer(4) oR ! recombination length option (1 = constant, 2 = exponential)
real(8) uRg ! recombination rate
real(8) uREg ! recombination error rate
real(8) lRm ! recombination length, mean
integer(4) lRmx ! recombination length, mean
integer(4) lR ! recombination length
integer(4) Rstart ! recombination start position
integer(4) Rend ! recombination start position
integer(4) krec ! recombination flag
real(8) fR ! factor on recombination rate (temporary)
integer(4) xR ! multiplier on recombinations
integer(4) xRuse

!
! -------------
! --- death ---
! -------------
!
real(8) kd ! death rate

!
! --------------
! --- output ---
! --------------
!
! --- general ---
!
character(4) sap

!
! --- NOTE ---
!
real(8) dtNote
real(8) tNoteNext

!
! --- STA ---
!
integer(4) opSTA ! output option (species calc: 0 = no, 1 = yes)
integer(4) oFixOut ! fixed changes output, 0 = no, 1 = yes
integer(4) opSTAf ! print STA at end option, 0 = no, 1 = yes
real(8) tpsSTA ! output start
real(8) dtpSTA ! output time step
real(8) tpnSTA ! next output time
integer(4) ipSTAf ! final print flag 
integer(4) ipSTAn ! counter
!
! total counters
!
real(8) nadt, nsrdt ! number of divisions
real(8), allocatable :: nadt_sp(:), nsrdt_sp(:)
real(8) fmutt, frect, fmutrecn ! fraction mut and reco at dilution
real(8), allocatable :: fmutt_sp(:)
real(8), allocatable :: frect_sp(:)
real(8), allocatable :: fmutrecn_sp(:)
real(8) namt, nsrmt ! number of mutations
real(8), allocatable :: namt_sp(:), nsrmt_sp(:)
real(8) nmTt, nmSt, nmNt, nmXt
real(8) nmTij(4,4), nmSij(4,4), nmNij(4,4), nmXij(4,4)
real(8), allocatable :: nmTt_sp(:)
real(8), allocatable :: nmSt_sp(:)
real(8), allocatable :: nmNt_sp(:)
real(8), allocatable :: nmXt_sp(:)
real(8), allocatable :: nmTij_sp(:,:,:)
real(8), allocatable :: nmSij_sp(:,:,:)
real(8), allocatable :: nmNij_sp(:,:,:)
real(8), allocatable :: nmXij_sp(:,:,:)
real(8) nart, nsrrt ! number of recombinations
real(8), allocatable :: nart_sp(:), nsrrt_sp(:)
real(8) nrTt, nrSt, nrNt, nrXt
real(8) nrTij(4,4), nrSij(4,4), nrNij(4,4), nrXij(4,4)
real(8), allocatable :: nrTt_sp(:)
real(8), allocatable :: nrSt_sp(:)
real(8), allocatable :: nrNt_sp(:)
real(8), allocatable :: nrXt_sp(:)
real(8), allocatable :: nrTij_sp(:,:,:)
real(8), allocatable :: nrSij_sp(:,:,:)
real(8), allocatable :: nrNij_sp(:,:,:)
real(8), allocatable :: nrXij_sp(:,:,:)
real(8) rdeltat, rnut, rn ! recombination statistics
real(8) rdelta, rnu
integer(4) irtype
real(8) rtype(10), rtypet
real(8), allocatable :: rdeltat_sp(:)
real(8), allocatable :: rnut_sp(:)
real(8), allocatable :: rn_sp(:)
real(8), allocatable :: rtype_sp(:,:)
real(8) fpt, fpn2, fp ! nonsynonymous penalty factor statistics
real(8), allocatable :: fpt_sp(:)
real(8), allocatable :: fpn2_sp(:)
real(8) nawt, nsrwt ! number of deaths/washouts
real(8), allocatable :: nawt_sp(:), nsrwt_sp(:)

!
! population averages etc
!
real(8) nda, cna, tga, kga
real(8) nma(4,4), nmTa, nmXa, nmSa, nmXSa, nmNa
real(8) nra(4,4), nrTa, nrXa, nrSa, nrXSa, nrNa
real(8) VCa, VNa, VPa
real(8) VmaxCa, VmaxNa, VmaxPa
real(8) nna(4)
real(8) fnta(4)
real(8) fntta(4)
real(8) faaa(20)
real(8) faastd(20)
real(8) naaxa(20,20)
real(8) qCa, qNa, qPa
real(8) q0Ca, q0Na, q0Pa
real(8) q0DNACa, q0DNANa, q0DNAPa
real(8) q0aaCa, q0aaNa
real(8) CLimita
real(8) NLimita
real(8) PLimita

real(8) kgave
real(8) kgmin
real(8) kgmax
real(8) tgstd
real(8) kgstd
!real(8) nmNave
integer(4) nmNmin
integer(4) nmNmax
real(8) nmNstd

integer(4) a_nmXt, a_nmXSt
real(8) a_GCTt, a_GCXt, a_GCSt, a_GCXSt, a_GCNt

real(8) CT_nmT_nmT
real(8) CT_nmT_nmX
real(8) CT_nmT_nmS
real(8) CT_nmT_nmXS
real(8) CT_nmT_nmN
real(8) CT_nmX_nmT
real(8) CT_nmX_nmX
real(8) CT_nmX_nmS
real(8) CT_nmX_nmXS
real(8) CT_nmX_nmN
real(8) CT_nmS_nmT
real(8) CT_nmS_nmX
real(8) CT_nmS_nmS
real(8) CT_nmS_nmXS
real(8) CT_nmS_nmN
real(8) CT_nmXS_nmT
real(8) CT_nmXS_nmX
real(8) CT_nmXS_nmS
real(8) CT_nmXS_nmXS
real(8) CT_nmXS_nmN
real(8) CT_nmN_nmT
real(8) CT_nmN_nmX
real(8) CT_nmN_nmS
real(8) CT_nmN_nmXS
real(8) CT_nmN_nmN
real(8) CT_nmT
real(8) CT_nmX
real(8) CT_nmS
real(8) CT_nmXS
real(8) CT_nmN

real(8) CT_GCT_GCT
real(8) CT_GCT_GCX
real(8) CT_GCT_GCS
real(8) CT_GCT_GCXS
real(8) CT_GCT_GCN
real(8) CT_GCX_GCT
real(8) CT_GCX_GCX
real(8) CT_GCX_GCS
real(8) CT_GCX_GCXS
real(8) CT_GCX_GCN
real(8) CT_GCS_GCT
real(8) CT_GCS_GCX
real(8) CT_GCS_GCS
real(8) CT_GCS_GCXS
real(8) CT_GCS_GCN
real(8) CT_GCXS_GCT
real(8) CT_GCXS_GCX
real(8) CT_GCXS_GCS
real(8) CT_GCXS_GCXS
real(8) CT_GCXS_GCN
real(8) CT_GCN_GCT
real(8) CT_GCN_GCX
real(8) CT_GCN_GCS
real(8) CT_GCN_GCXS
real(8) CT_GCN_GCN
real(8) CT_GCT
real(8) CT_GCX
real(8) CT_GCS
real(8) CT_GCXS
real(8) CT_GCN

real(8) AT_nmT
real(8) AT_nmX
real(8) AT_nmS
real(8) AT_nmXS
real(8) AT_nmN

real(8) AT_GCT
real(8) AT_GCX
real(8) AT_GCS
real(8) AT_GCXS
real(8) AT_GCN


!
! GCT, GCS, GCN, GCX, GCXS
!
real(8) GCT
real(8) GCS
real(8) GCN
real(8) GCX
real(8) GCXS
real(8) ntnka(4)
real(8) ntnSka(4)
real(8) ntnNka(4)
real(8) ntnXka(4)
real(8) ntnXkax(4)

!
! strain summary
!
integer(4), allocatable :: srsn(:), sn(:)
integer(4) is1, ns

!
! change statistics
!
real(8) t_first, t_last, t_now
real(8) GC_first, GC_last, GC_now
real(8) GCT_first, GCT_last, GCT_now
real(8) GCS_first, GCS_last, GCS_now
real(8) GCN_first, GCN_last, GCN_now
real(8) GCX_first, GCX_last, GCX_now
real(8) GCXS_first, GCXS_last, GCXS_now
real(8) piS_first, piS_last, piS_now
real(8) piN_first, piN_last, piN_now
real(8) piM_first, piM_last, piM_now
real(8) piR_first, piR_last, piR_now
real(8) piT_first, piT_last, piT_now
real(8) dNdS_first, dNdS_last, dNdS_now
real(8) rm_first, rm_last, rm_now
real(8) q0DNAC_first, q0DNAC_last, q0DNAC_now
real(8) q0DNAN_first, q0DNAN_last, q0DNAN_now
real(8) q0DNAP_first, q0DNAP_last, q0DNAP_now
real(8) q0aaC_first, q0aaC_last, q0aaC_now
real(8) q0aaN_first, q0aaN_last, q0aaN_now
real(8) q0C_first, q0C_last, q0C_now
real(8) q0N_first, q0N_last, q0N_now
real(8) q0P_first, q0P_last, q0P_now
real(8) kg_first, kg_last, kg_now

!
! --- GLO ---
!
integer(4) opGLO ! GLO output option: 0 = no, 1 = random npGLO, 2 = first npGLO
real(8) tpsGLO ! GLO start time
real(8) dtpGLO ! GLO output step
integer(4) npGLO ! no. cells to output 
real(8) tpnGLO ! next GLO output time
integer(4) ipGLO ! GLO output counter

!
! --- CLK ---
!
integer(4) opCLK ! output CLK? 0 = no, 1 = yes
integer(4) ipCLK ! CLK output flag

!
! --- POP ---
!
integer(4) opPOP ! POP output option: 0 = none, 1 = specified steps - not used
real(8) tpsPOP ! POP start time
real(8) dtpPOP ! POP output step
real(8) tpnPOP ! next POP output time
integer(4) ipPOP ! POP output counter

!
! --- ORF ---
!
integer(4) opORF ! ORF output option: -1 = none, 0 = all, 1 = random npORF, 2 = first npORF
real(8) tpsORF ! ORF start time
real(8) dtpORF ! ORF output step
integer(4) npORF ! no. cells to output
real(8) tpnORF ! next ORF output time

real(8) SRTp
real(8) pnmta, pnmSa, pnmNa
real(8) pnrta, pnrSa, pnrNa
real(8) ppiSa, ppiS
real(8) ppiNa, ppiN
real(8) ppiMa, ppiM 
real(8) ppiRa, ppiR 
real(8) ppiTa, ppiT

!
! ---------------------------
! --- time-variable input ---
! ---------------------------
!
integer(4) it ! time-variable input index
real(8) fQ
real(8) fSinC
real(8) fSinN
real(8) fSinP
real(8) fSCfix
real(8) fSNfix
real(8) fSPfix
real(8) fX1
real(8) fX2
real(8), dimension(:), allocatable :: tt ! time-variable input times
real(8), dimension(:), allocatable :: tfQ
real(8), dimension(:), allocatable :: tfSinC
real(8), dimension(:), allocatable :: tfSinN
real(8), dimension(:), allocatable :: tfSinP
real(8), dimension(:), allocatable :: tfSCfix
real(8), dimension(:), allocatable :: tfSNfix
real(8), dimension(:), allocatable :: tfSPfix
real(8), dimension(:), allocatable :: tfX1
real(8), dimension(:), allocatable :: tfX2

!
! -----------------
! -----------------
! --- start msg ---
! -----------------
! -----------------
!
write(*,*) 'IAM:'
write(*,*) 'IAM: IAM START'
write(*,*) 'IAM:'

!
! ---------------
! ---------------
! --- formats ---
! ---------------
! ---------------
!
!         1         2         3         4         5
!12345678901234567890123456789012345678901234567890123456789
!>IAMGLO00010001 IAM GENOME GLO ip: 0001 ia: 0001.
205 format(1a49) ! DNA header, write
201 format(80a1) ! DNA, write
206 format(1e24.16,7i10,12i10,1e24.16,3i10) ! VAR, write
208 format(1a100) ! DNA header, read
207 format(70a1) ! DNA, read

!
! ------------------------
! ------------------------
! --- read const input ---
! ------------------------
! ------------------------
!
! -----------------
! --- open file ---
! -----------------
!
open(111,file='IAM_I_c.txt')

!
! ------------------
! --- dimensions ---
! ------------------
!
read(111,*) dimna
write(*,*) 'IAM\INPUT: dimna = ', dimna 
read(111,*) dimnc
write(*,*) 'IAM\INPUT: dimnc = ', dimnc
read(111,*) dimng
write(*,*) 'IAM\INPUT: dimng  = ', dimng
read(111,*) dimnnt
write(*,*) 'IAM\INPUT: dimnnt = ', dimnnt
read(111,*) dimnp
write(*,*) 'IAM\INPUT: dimnp = ', dimnp
read(111,*) dimnsp
write(*,*) 'IAM\INPUT: dimnsp = ', dimnsp
read(111,*) dimntb
write(*,*) 'IAM\INPUT: dimntb = ', dimntb
read(111,*) dimns
write(*,*) 'IAM\INPUT: dimns = ', dimns
read(111,*) dimnk
write(*,*) 'IAM\INPUT: dimnk = ', dimnk
read(111,*) ! SPARE
read(111,*) cTime
write(*,*) 'IAM\INPUT: cTime = ', cTime
read(111,*) cComp
write(*,*) 'IAM\INPUT: cComp = ', cComp
read(111,*) cIter
write(*,*) 'IAM\INPUT: cIter = ', cIter
if (cIter.lt.1) then ! legacy support
    cIter = 1
endif
read(111,*) ! SPARE

!
! -------------------
! --- run control ---
! -------------------
!
read(111,*) tend
write(*,*) 'IAM\INPUT: tend = ', tend
read(111,*) dt
write(*,*) 'IAM\INPUT: dt = ', dt
read(111,*) oHot
write(*,*) 'IAM\INPUT: oHot = ', oHot
read(111,*) oHANS
write(*,*) 'IAM\INPUT: oHANS = ', oHANS
read(111,*) oHANSu
write(*,*) 'IAM\INPUT: oHANSu = ', oHANSu
read(111,*) oHANSd
write(*,*) 'IAM\INPUT: oHANSd = ', oHANSd
read(111,*) oHANSs
write(*,*) 'IAM\INPUT: oHANSs = ', oHANSs
read(111,*) RSi
write(*,*) 'IAM\INPUT: RSi = ', RSi
read(111,*) RSj
write(*,*) 'IAM\INPUT: RSj = ', RSj
read(111,*) oHots
if (oHots.lt.0.) then ! legacy support
    oHots = 1
endif
write(*,*) 'IAM\INPUT: oHots = ', oHots
read(111,*) oHANSr
write(*,*) 'IAM\INPUT: oHANSr = ', oHANSr
if (oHANSr.lt.0.) then ! legacy support
    oHANSr = 1
endif
read(111,*) oIter
write(*,*) 'IAM\INPUT: oIter = ', oIter
read(111,*) nIter
write(*,*) 'IAM\INPUT: nIter = ', nIter
read(111,*) ! SPARE
!write(*,*) 'IAM\INPUT: SPARE = ', SPARE

!
! -------------
! --- cells ---
! -------------
!
read(111,*) oSlim
write(*,*) 'IAM\INPUT: oSlim = ', oSlim
read(111,*) Popic
write(*,*) 'IAM\INPUT: Popic = ', Popic
read(111,*) ! SPARE
!write(*,*) 'IAM\INPUT: SPARE = ', SPARE
read(111,*) ! SPARE
!write(*,*) 'IAM\INPUT: SPARE = ', SPARE
read(111,*) ! SPARE
!write(*,*) 'IAM\INPUT: SPARE = ', SPARE
read(111,*) ! SPARE
!write(*,*) 'IAM\INPUT: SPARE = ', SPARE
read(111,*) ! SPARE
!write(*,*) 'IAM\INPUT: SPARE = ', SPARE

!
! -------------------------
! --- genome & proteins ---
! -------------------------
!
read(111,*) oStartDNA
write(*,*) 'IAM\INPUT: oStartDNA = ', oStartDNA
read(111,*) ipw
write(*,*) 'IAM\INPUT: ipw = ', ipw
read(111,*) nger
write(*,*) 'IAM\INPUT: nger = ', nger
read(111,*) ficA
write(*,*) 'IAM\INPUT: ficA = ', ficA
read(111,*) ficT
write(*,*) 'IAM\INPUT: ficT = ', ficT
read(111,*) ficC
write(*,*) 'IAM\INPUT: ficC = ', ficC
read(111,*) ficG
write(*,*) 'IAM\INPUT: ficG = ', ficG
read(111,*) pminspacing
write(*,*) 'IAM\INPUT: pminspacing = ', pminspacing
read(111,*) pmaxspacing
write(*,*) 'IAM\INPUT: pmaxspacing = ', pmaxspacing
read(111,*) pminlength
write(*,*) 'IAM\INPUT: pminlength = ', pminlength
read(111,*) pmaxlength
write(*,*) 'IAM\INPUT: pmaxlength = ', pmaxlength
read(111,*) pfAT
write(*,*) 'IAM\INPUT: pfAT = ', pfAT
if (pfAT.lt.0.) then ! legacy support
    pfAT = 0.
endif
read(111,*) pfGC
write(*,*) 'IAM\INPUT: pfGC = ', pfGC
if (pfGC.lt.0.) then ! legacy support
    pfGC = 0.
endif
read(111,*) ! SPARE
!write(*,*) 'IAM\INPUT: SPARE = ', SPARE
read(111,*) ! SPARE
!write(*,*) 'IAM\INPUT: SPARE = ', SPARE
read(111,*) ! SPARE
!write(*,*) 'IAM\INPUT: SPARE = ', SPARE

!
! --------------------------
! --- reactor simulation ---
! --------------------------
!
read(111,*) oFixExt
write(*,*) 'IAM\INPUT: oFixExt = ', oFixExt
if (oHANS.eq.1) then
    if (oFixExt.eq.0) then
        write(*,*) 'IAM\HANS: Warning: oFixExt needs to be 1+.'
        write(*,*) 'IAM\HANS: Adjusting...'
        oFixExt = 1
    endif
endif
read(111,*) oImpExt
write(*,*) 'IAM\INPUT: oImpExt = ', oImpExt
read(111,*) SC0
write(*,*) 'IAM\INPUT: SC0 = ', SC0
read(111,*) SN0
write(*,*) 'IAM\INPUT: SN0 = ', SN0
read(111,*) SP0
write(*,*) 'IAM\INPUT: SP0 = ', SP0
read(111,*) SinCB
write(*,*) 'IAM\INPUT: SinCB = ', SinCB
read(111,*) SinNB
write(*,*) 'IAM\INPUT: SinNB = ', SinNB
read(111,*) SinPB
write(*,*) 'IAM\INPUT: SinPB = ', SinPB
SinC = SinCB
SinN = SinNB
SinP = SinPB
read(111,*) VOL
write(*,*) 'IAM\INPUT: VOL = ', VOL
VOL = VOL * 1e-3
read(111,*) QB
write(*,*) 'IAM\INPUT: QB = ', QB
QB = QB * 1e-3
Q = QB
read(111,*) oPD
write(*,*) 'IAM\INPUT: oPD = ', oPD
read(111,*) dtD
write(*,*) 'IAM\INPUT: dtD = ', dtD
read(111,*) fD
write(*,*) 'IAM\INPUT: fD = ', fD
read(111,*) oMix
write(*,*) 'IAM\INPUT: oMix = ', oMix
read(111,*) dtmix
write(*,*) 'IAM\INPUT: dtmix = ', dtmix
read(111,*) nspmix
write(*,*) 'IAM\INPUT: nspmix = ', nspmix
if (nspmix.le.0) then
    nspmix = dimnsp
endif
read(111,*) ! SPARE
!write(*,*) 'IAM\INPUT: SPARE = ', SPARE
read(111,*) ! SPARE
!write(*,*) 'IAM\INPUT: SPARE = ', SPARE
read(111,*) ! SPARE
!write(*,*) 'IAM\INPUT: SPARE = ', SPARE
read(111,*) ! SPARE
!write(*,*) 'IAM\INPUT: SPARE = ', SPARE

!
! ------------------------------------------------------
! --- nutrient uptake and intracellular mass balance ---
! ------------------------------------------------------
!
read(111,*) Vmax0C
write(*,*) 'IAM\INPUT: Vmax0C = ', Vmax0C
read(111,*) Vmax0N
write(*,*) 'IAM\INPUT: Vmax0N = ', Vmax0N
read(111,*) Vmax0P
write(*,*) 'IAM\INPUT: Vmax0P = ', Vmax0P
read(111,*) oLN
write(*,*) 'IAM\INPUT: oLN = ', oLN
read(111,*) KmC
write(*,*) 'IAM\INPUT: KmC = ', KmC
read(111,*) KmN
write(*,*) 'IAM\INPUT: KmN = ', KmN
read(111,*) KmP
write(*,*) 'IAM\INPUT: KmP = ', KmP
read(111,*) oLP
write(*,*) 'IAM\INPUT: oLP = ', oLP
read(111,*) K
write(*,*) 'IAM\INPUT: K = ', K
read(111,*) oLM
write(*,*) 'IAM\INPUT: oLM = ', oLM
read(111,*) oKimura
write(*,*) 'IAM\INPUT: oKimura = ', oKimura
read(111,*) fpn
write(*,*) 'IAM\INPUT: fpn = ', fpn
read(111,*) fps
write(*,*) 'IAM\INPUT: fps = ', fps
read(111,*) fpx
write(*,*) 'IAM\INPUT: fpx = ', fpx
read(111,*) rnn
write(*,*) 'IAM\INPUT: rnn = ', rnn
read(111,*) rns
write(*,*) 'IAM\INPUT: rns = ', rns
read(111,*) rnx
write(*,*) 'IAM\INPUT: rnx = ', rnx
read(111,*) fpbeta
write(*,*) 'IAM\INPUT: fpbeta = ', fpbeta
Mdc = 3.0 ! hardwired, make input
Mfc = rnn ! hardwired, make input
read(111,*) oLoth
write(*,*) 'IAM\INPUT: oLoth = ', oLoth
read(111,*) oLothIC
write(*,*) 'IAM\INPUT: oLothIC = ', oLothIC
read(111,*) sGCt
write(*,*) 'IAM\INPUT: sGCt = ', sGCt
read(111,*) fGCt0
write(*,*) 'IAM\INPUT: fGCt0 = ', fGCt0
read(111,*) oCAD
write(*,*) 'IAM\INPUT: oCAD = ', oCAD
read(111,*) oFTSS
write(*,*) 'IAM\INPUT: oFTSS = ', oFTSS
if (oLM.ge.5) then ! hardwired based on oLM
    oFTSS = 1
else
    oFTSS = 0
endif
!fdel = 0.99 ! hardwired, make input
read(111,*) fdel
write(*,*) 'IAM\INPUT: fdel = ', fdel
read(111,*) ! SPARE
!write(*,*) 'IAM\INPUT: SPARE = ', SPARE
read(111,*) ! SPARE
!write(*,*) 'IAM\INPUT: SPARE = ', SPARE

!
! ----------------
! --- division ---
! ----------------
!
read(111,*) fdq0
write(*,*) 'IAM\INPUT: fdq0 = ', fdq0
if (fdq0.lt.0.) then ! legacy support
    fdq0 = 1.
endif
read(111,*) oFixCells
write(*,*) 'IAM\INPUT: oFixCells = ', oFixCells
read(111,*) oLQ
write(*,*) 'IAM\INPUT: oLQ = ', oLQ
if (oLQ.eq.1) then ! legacy support
    oLQ = 31
endif
read(111,*) oKillDaughter
write(*,*) 'IAM\INPUT: oKillDaughter = ', oKillDaughter
read(111,*) RCVmm
write(*,*) 'IAM\INPUT: RCVmm = ', RCVmm
read(111,*) RCVmcv
write(*,*) 'IAM\INPUT: RCVmcv = ', RCVmcv
read(111,*) RCVmcvl
write(*,*) 'IAM\INPUT: RCVmcvl = ', RCVmcvl
read(111,*) RCVmll
write(*,*) 'IAM\INPUT: RCVmll = ', RCVmll
read(111,*) RCVmhl
write(*,*) 'IAM\INPUT: RCVmhl = ', RCVmhl
read(111,*) n0DNA
write(*,*) 'IAM\INPUT: n0DNA = ', n0DNA
read(111,*) m0aa
write(*,*) 'IAM\INPUT: m0aa = ', m0aa
read(111,*) q0othC
write(*,*) 'IAM\INPUT: q0othC = ', q0othC
read(111,*) q0othN
write(*,*) 'IAM\INPUT: q0othN = ', q0othN
read(111,*) q0othP
write(*,*) 'IAM\INPUT: q0othP = ', q0othP
read(111,*) fqmaxC
write(*,*) 'IAM\INPUT: fqmaxC = ', fqmaxC
read(111,*) fqmaxN
write(*,*) 'IAM\INPUT: fqmaxN = ', fqmaxN
read(111,*) fqmaxP
write(*,*) 'IAM\INPUT: fqmaxP = ', fqmaxP
read(111,*) ! SPARE
!write(*,*) 'IAM\INPUT: SPARE = ', SPARE
read(111,*) ! SPARE
!write(*,*) 'IAM\INPUT: SPARE = ', SPARE
read(111,*) ! SPARE
!write(*,*) 'IAM\INPUT: SPARE = ', SPARE
read(111,*) ! SPARE
!write(*,*) 'IAM\INPUT: SPARE = ', SPARE
read(111,*) ! SPARE
!write(*,*) 'IAM\INPUT: SPARE = ', SPARE

!
! --------------------------------
! --- mutation & recombination ---
! --------------------------------
!
read(111,*) fu
write(*,*) 'IAM\INPUT: fu = ', fu
ubij = 0.
read(111,*) ubij(1,2)
read(111,*) ubij(1,3)
read(111,*) ubij(1,4)
read(111,*) ubij(2,1)
read(111,*) ubij(2,3)
read(111,*) ubij(2,4)
read(111,*) ubij(3,1)
read(111,*) ubij(3,2)
read(111,*) ubij(3,4)
read(111,*) ubij(4,1)
read(111,*) ubij(4,2)
read(111,*) ubij(4,3)
do nt1 = 1, 4
    do nt2 = 1, 4
        if (nt1.ne.nt2) then
            write(*,*) 'IAM\INPUT: ubij, i, j = ', ubij(nt1,nt2), nt1, nt2
        endif
    enddo
enddo
read(111,*) fM
write(*,*) 'IAM\INPUT: fM = ', fM
read(111,*) oR 
write(*,*) 'IAM\INPUT: oR = ', oR
read(111,*) uRg
write(*,*) 'IAM\INPUT: uRg = ', uRg
read(111,*) lRm
write(*,*) 'IAM\INPUT: lRm = ', lRm
read(111,*) lRmx
write(*,*) 'IAM\INPUT: lRmx = ', lRmx
read(111,*) fR
write(*,*) 'IAM\INPUT: fR = ', fR
read(111,*) xM
write(*,*) 'IAM\INPUT: xM = ', xM
read(111,*) xR
write(*,*) 'IAM\INPUT: xR = ', xR
read(111,*) fmKill
write(*,*) 'IAM\INPUT: fmKill = ', fmKill
read(111,*) ! SPARE
!write(*,*) 'IAM\INPUT: SPARE = ', SPARE
read(111,*) ! SPARE
!write(*,*) 'IAM\INPUT: SPARE = ', SPARE

!
! -------------
! --- death ---
! -------------
!
read(111,*) kd
write(*,*) 'IAM\INPUT: kd = ', kd
read(111,*) ! SPARE
!write(*,*) 'IAM\INPUT: SPARE = ', SPARE
read(111,*) ! SPARE
!write(*,*) 'IAM\INPUT: SPARE = ', SPARE
read(111,*) ! SPARE
!write(*,*) 'IAM\INPUT: SPARE = ', SPARE
read(111,*) ! SPARE
!write(*,*) 'IAM\INPUT: SPARE = ', SPARE
read(111,*) ! SPARE
!write(*,*) 'IAM\INPUT: SPARE = ', SPARE

!
! --------------
! --- output ---
! --------------
!
read(111,*) dtNote
write(*,*) 'IAM\INPUT: dtNote = ', dtNote
read(111,*) opSTA
write(*,*) 'IAM\INPUT: opSTA = ', opSTA
read(111,*) oFixOut
write(*,*) 'IAM\INPUT: oFixOut = ', oFixOut
opSTAf = 1 ! make input
read(111,*) tpsSTA
write(*,*) 'IAM\INPUT: tpsSTA = ', tpsSTA
read(111,*) dtpSTA
write(*,*) 'IAM\INPUT: dtpSTA = ', dtpSTA
read(111,*) opGLO
write(*,*) 'IAM\INPUT: opGLO = ', opGLO
read(111,*) tpsGLO
write(*,*) 'IAM\INPUT: tpsGLO = ', tpsGLO
read(111,*) dtpGLO
write(*,*) 'IAM\INPUT: dtpGLO = ', dtpGLO
read(111,*) npGLO
write(*,*) 'IAM\INPUT: npGLO = ', npGLO
read(111,*) opCLK
write(*,*) 'IAM\INPUT: opCLK = ', opCLK
read(111,*) opPOP
write(*,*) 'IAM\INPUT: opPOP = ', opPOP
read(111,*) tpsPOP
write(*,*) 'IAM\INPUT: tpsPOP = ', tpsPOP
read(111,*) dtpPOP
write(*,*) 'IAM\INPUT: dtpPOP = ', dtpPOP
read(111,*) opORF
write(*,*) 'IAM\INPUT: opORF = ', opORF
read(111,*) tpsORF
write(*,*) 'IAM\INPUT: tpsORF = ', tpsORF
read(111,*) dtpORF
write(*,*) 'IAM\INPUT: dtpORF = ', dtpORF
read(111,*) npORF
write(*,*) 'IAM\INPUT: npORF = ', npORF
read(111,*) !SPARE
!write(*,*) 'IAM\INPUT: SPARE = ', SPARE

!
! ------------------
! --- close file ---
! ------------------
!
close(111)

!
! ---------------------------
! ---------------------------
! --- parallel processing ---
! ---------------------------
! ---------------------------
!
ntd = dimnsp
write(*,*) 'IAM\PAR: Using ntd = ', ntd

!
! --------------------------
! --------------------------
! --- allocate variables ---
! --------------------------
! --------------------------
!
write(*,*) 'IAM\PAR: Allocating variables. Start.'

dimna_sp = dble(dimna) / dble(dimnsp)

allocate (a_sa(2,dimnsp,dimna_sp))
allocate (a_sr(2,dimnsp,dimna_sp))
allocate (a_srm(2,dimnsp,dimna_sp))
allocate (a_srr(2,dimnsp,dimna_sp))
allocate (a_iar(2,dimnsp,dimna_sp))
allocate (a_gid(2,dimnsp,dimna_sp))
allocate (a_cn(2,dimnsp,dimna_sp))
allocate (a_ck(dimnc,dimnsp,dimna_sp,2))
allocate (a_cnt(dimnc,dimnsp,dimna_sp,2))
allocate (a_cfp(dimnc,dimnsp,dimna_sp,2))
allocate (a_ntnk(2,dimnsp,dimna_sp,4))
allocate (a_naa(2,dimnsp,dimna_sp,20))
allocate (a_naax(2,dimnsp,dimna_sp,20,20))
allocate (a_uEg(2,dimnsp,dimna_sp))
allocate (a_q0DNAC(2,dimnsp,dimna_sp))
allocate (a_q0DNAN(2,dimnsp,dimna_sp))
allocate (a_q0DNAP(2,dimnsp,dimna_sp))
allocate (a_q0aaC(2,dimnsp,dimna_sp))
allocate (a_q0aaN(2,dimnsp,dimna_sp))
allocate (a_VmaxC(2,dimnsp,dimna_sp))
allocate (a_VmaxN(2,dimnsp,dimna_sp))
allocate (a_VmaxP(2,dimnsp,dimna_sp))
if (oSlim.eq.0) then
    allocate (a_sn(2,dimnsp,dimna_sp))
    allocate (a_rsn(2,dimnsp,dimna_sp))
    allocate (a_nd(2,dimnsp,dimna_sp))
    allocate (a_ntnSk(2,dimnsp,dimna_sp,4))
    allocate (a_ntnNk(2,dimnsp,dimna_sp,4))
    allocate (a_ntnXk(2,dimnsp,dimna_sp,4))
    allocate (a_nmt(2,dimnsp,dimna_sp))
    allocate (a_nmS(2,dimnsp,dimna_sp))
    allocate (a_nmN(2,dimnsp,dimna_sp))
    allocate (a_nrt(2,dimnsp,dimna_sp))
    allocate (a_nrS(2,dimnsp,dimna_sp))
    allocate (a_nrN(2,dimnsp,dimna_sp))
    allocate (a_cXSN(dimnc,dimnsp,dimna_sp,2))
    allocate (a_cmr(dimnc,dimnsp,dimna_sp,2))
    allocate (a_nm(2,dimnsp,dimna_sp,4,4))
    allocate (a_nr(2,dimnsp,dimna_sp,4,4))
    allocate (a_nntt(2,dimnsp,dimna_sp,4))
    allocate (a_VC(2,dimnsp,dimna_sp))
    allocate (a_VN(2,dimnsp,dimna_sp))
    allocate (a_VP(2,dimnsp,dimna_sp))
    allocate (a_qC(2,dimnsp,dimna_sp))
    allocate (a_qN(2,dimnsp,dimna_sp))
    allocate (a_qP(2,dimnsp,dimna_sp))
    allocate (a_CLimit(2,dimnsp,dimna_sp))
    allocate (a_NLimit(2,dimnsp,dimna_sp))
    allocate (a_PLimit(2,dimnsp,dimna_sp))
    allocate (a_tb(2,dimnsp,dimna_sp))
    allocate (a_tg(2,dimnsp,dimna_sp))
endif

allocate (g_nt(dimng,dimnnt))
allocate (g_nt_o(dimng,dimnnt))
allocate (g_ntfS(dimng,dimnnt))
allocate (g_ntn(dimng))
allocate (g_ntnS(dimng))
allocate (g_ntnN(dimng))
allocate (g_ntnX(dimng))
allocate (g_ntnk(dimng,4))
allocate (g_ntnSk(dimng,4))
allocate (g_ntnNk(dimng,4))
allocate (g_ntnXk(dimng,4))
allocate (g_ntnij(dimng,4,4))
allocate (g_ntnSij(dimng,4,4))
allocate (g_ntnNij(dimng,4,4))
allocate (g_ntnXij(dimng,4,4))
allocate (g_nntt(dimng,4))
allocate (g_nnttt(dimng))
allocate (g_nnttw(dimng,4))
allocate (g_nntttw(dimng))
allocate (g_naa(dimng,20))
allocate (g_uEg(dimng))
allocate (g_ip(dimng,dimnnt))
allocate (g_ntu1(dimnsp,dimnnt))
allocate (g_ntu2(dimnsp,dimnnt))
allocate (g_cmru1(dimnsp,dimnnt))
allocate (g_cmru2(dimnsp,dimnnt))
allocate (g_cic1(dimnsp,dimnnt))
!allocate (g_cic2(dimnsp,dimnnt))
allocate (g_cfp1(dimnsp,dimnnt))
allocate (g_cfp2(dimnsp,dimnnt))
allocate (g_ncf(dimng))
allocate (g_ncSf(dimng))
allocate (g_ncNf(dimng))
allocate (g_ncmf(dimng))
allocate (g_ncrf(dimng))
allocate (g_ncff(dimng))
allocate (g_ncSff(dimng))
allocate (g_ncNff(dimng))
allocate (g_ncmff(dimng))
allocate (g_ncrff(dimng))
allocate (g_cmr(dimng,dimnnt))
allocate (g_q0DNAC(dimng))
allocate (g_q0DNAN(dimng))
allocate (g_q0DNAP(dimng))
allocate (g_q0aaC(dimng))
allocate (g_q0aaN(dimng))
allocate (g_Vmax0C(dimng))
allocate (g_Vmax0N(dimng))
allocate (g_Vmax0P(dimng))
allocate (g_pntstart(dimng,dimnp))
allocate (g_pntstop(dimng,dimnp))
allocate (g_pother(dimng,dimnp))
allocate (g_pgamma(dimng,dimnp))
allocate (g_pntn(dimng,dimnp))
allocate (g_pntnS(dimng,dimnp))
allocate (g_pntnN(dimng,dimnp))
allocate (g_pnmS(dimng,dimnp))
allocate (g_pnmN(dimng,dimnp))
allocate (g_pnrS(dimng,dimnp))
allocate (g_pnrN(dimng,dimnp))
allocate (g_np(dimng))
allocate (g_naat(dimng))

allocate (ifx(dimng))
allocate (fx_cn(dimng))
allocate (fx_ck(dimnc,dimng))
allocate (fx_cnt(dimnc,dimng))
allocate (fx_cXSN(dimnc,dimng))
allocate (fx_cmr(dimnc,dimng))
allocate (ifx_sp(dimnsp,dimng))
allocate (fx_cn_sp(dimnsp,dimng))
allocate (fx_ck_sp(dimnc,dimnsp,dimng))
allocate (fx_cnt_sp(dimnc,dimnsp,dimng))
allocate (fx_cXSN_sp(dimnc,dimnsp,dimng))
allocate (fx_cmr_sp(dimnc,dimnsp,dimng))

allocate (t_sp(dimnsp))
allocate (nat_comp_sp(100,dimnsp))
allocate (SRT_comp_sp(100,dimnsp))
allocate (err_sp(dimnsp))
allocate (idum_sp(dimnsp))
allocate (cp_ENS_t1_sp(dimnsp), cp_ENS_delta_sp(dimnsp))
allocate (cp_HANS1_t1_sp(dimnsp), cp_HANS1_delta_sp(dimnsp))
allocate (cp_HANS2_t1_sp(dimnsp), cp_HANS2_delta_sp(dimnsp))
!allocate (cp_HANS2a_t1_sp(dimnsp), cp_HANS2a_delta_sp(dimnsp))
!allocate (cp_HANS2b_t1_sp(dimnsp), cp_HANS2b_delta_sp(dimnsp))
!allocate (cp_HANS2c_t1_sp(dimnsp), cp_HANS2c_delta_sp(dimnsp))
allocate (cp_mr1_t1_sp(dimnsp), cp_mr1_delta_sp(dimnsp))
allocate (cp_mr2_t1_sp(dimnsp), cp_mr2_delta_sp(dimnsp))
!allocate (cp_mr2a_t1_sp(dimnsp), cp_mr2a_delta_sp(dimnsp))
!allocate (cp_mr2b_t1_sp(dimnsp), cp_mr2b_delta_sp(dimnsp))
!allocate (cp_mr2c_t1_sp(dimnsp), cp_mr2c_delta_sp(dimnsp))
allocate (nat_sp(2,dimnsp))
allocate (naa_sp(dimnsp))
allocate (naf_sp(dimnsp))
allocate (iaf_sp(dimnsp,dimna_sp))
allocate (nar_sp(dimnsp))
allocate (iad_sp(dimnsp,dimna_sp*xR))
allocate (SRT_sp(dimnsp))
allocate (gsn_sp(dimnsp))
allocate (Popic_sp(dimnsp))
allocate (ncg_sp(dimnsp))
allocate (SC_sp(dimnsp))
allocate (SN_sp(dimnsp))
allocate (SP_sp(dimnsp))
allocate (RSC_sp(dimnsp))
allocate (RSN_sp(dimnsp))
allocate (RSP_sp(dimnsp))
allocate (VOL_sp(dimnsp))
allocate (Q_sp(dimnsp))
allocate (tDNext_sp(dimnsp))
allocate (VmaxCa_sp(dimnsp))
allocate (VmaxNa_sp(dimnsp))
allocate (VmaxPa_sp(dimnsp))
allocate (nadt_sp(dimnsp))
allocate (nsrdt_sp(dimnsp))
allocate (fmutt_sp(dimnsp))
allocate (frect_sp(dimnsp))
allocate (fmutrecn_sp(dimnsp))
allocate (namt_sp(dimnsp))
allocate (nsrmt_sp(dimnsp))
allocate (nmTt_sp(dimnsp))
allocate (nmSt_sp(dimnsp))
allocate (nmNt_sp(dimnsp))
allocate (nmXt_sp(dimnsp))
allocate (nmTij_sp(dimnsp,4,4))
allocate (nmSij_sp(dimnsp,4,4))
allocate (nmNij_sp(dimnsp,4,4))
allocate (nmXij_sp(dimnsp,4,4))
allocate (nart_sp(dimnsp))
allocate (nsrrt_sp(dimnsp))
allocate (nrTt_sp(dimnsp))
allocate (nrSt_sp(dimnsp))
allocate (nrNt_sp(dimnsp))
allocate (nrXt_sp(dimnsp))
allocate (nrTij_sp(dimnsp,4,4))
allocate (nrSij_sp(dimnsp,4,4))
allocate (nrNij_sp(dimnsp,4,4))
allocate (nrXij_sp(dimnsp,4,4))
allocate (rdeltat_sp(dimnsp))
allocate (rnut_sp(dimnsp))
allocate (rn_sp(dimnsp))
allocate (rtype_sp(dimnsp,10))
allocate (fpt_sp(dimnsp))
allocate (fpn2_sp(dimnsp))
allocate (nawt_sp(dimnsp))
allocate (nsrwt_sp(dimnsp))

allocate (tt(dimntb))
allocate (tfQ(dimntb))
allocate (tfSinC(dimntb))
allocate (tfSinN(dimntb))
allocate (tfSinP(dimntb))
allocate (tfSCfix(dimntb))
allocate (tfSNfix(dimntb))
allocate (tfSPfix(dimntb))
allocate (tfX1(dimntb))
allocate (tfX2(dimntb))

allocate (sntr(70))

allocate (sarkarm(dimns))
allocate (sarkarr(dimns))

allocate (kimura(dimnk))

allocate (srsn(dimna))
allocate (sn(dimna))

write(*,*) 'IAM\PAR: Allocating variables. End.'

!
! ----------------------------------
! ----------------------------------
! --- random number seed & setup ---
! ----------------------------------
! ----------------------------------
!
if (RSi.eq.-1.) then
    idum = -secnds(0d0)
else
    idum = -abs(RSi)
end if
if (cIter.gt.1) then
    idum = idum + cIter * dimnsp
endif
r = r4_uni(idum)
write(*,*) 'IAM\RAN: Random number seed i = ', idum, r
if (RSj.eq.-1.) then
    jdum = -secnds(0d0)
else
    jdum = -abs(RSj)
end if
!if (cIter.gt.1) then
!    jdum = jdum + cIter * dimnsp
!endif
r = r4_uni(jdum)
write(*,*) 'IAM\RAN: Random number seed j = ', jdum, r
do isp1 = 1, dimnsp
    idum_sp(isp1) = idum - isp1
enddo

call r4_nor_setup(kn,fn,wn)

!rfirst = .true.
!r = random_gamma(0.5,rfirst)
!rfirst = .false.

!
! ---------------------
! --- kimura preset ---
! ---------------------
!
!if ((oLM.eq.3).and.(oKimura.eq.2)) then
!    write(*,*) 'IAM\KIMURA: Populating preset arrays.'
!    do ia1 = 1, dimnk
!        r = random_gamma(fpbeta,rfirst)
!        kimura(ia1) = r
!    enddo        
!    write(*,*) 'IAM\KIMURA: Done.'
!endif


!
! --------------------
! --------------------
! --- prelim calcs ---
! --------------------
! --------------------
!
!
int_opo(1) = 2
int_opo(2) = 1
int_opo(3) = 4
int_opo(4) = 3

!
! -----------------
! --- nt lookup ---
! -----------------
!
int_snt(1) = 'A'
int_snt(2) = 'T'
int_snt(3) = 'C'
int_snt(4) = 'G'

cint(1,1) = 2
cint(1,2) = 3
cint(1,3) = 4
cint(2,1) = 1
cint(2,2) = 3
cint(2,3) = 4
cint(3,1) = 1
cint(3,2) = 2
cint(3,3) = 4
cint(4,1) = 1
cint(4,2) = 2
cint(4,3) = 3

!
! --------------------
! --- genetic code ---
! --------------------
!
gc = -1

gc(124) = 21
    
gc(211) = 22
gc(241) = 22
gc(214) = 22

gc(432) = 1
gc(433) = 1
gc(431) = 1
gc(434) = 1
    
gc(342) = 2
gc(343) = 2
gc(341) = 2
gc(344) = 2
gc(141) = 2
gc(144) = 2
    
gc(112) = 3
gc(113) = 3
    
gc(412) = 4
gc(413) = 4
    
gc(242) = 5
gc(243) = 5
    
gc(311) = 6
gc(314) = 6
    
gc(411) = 7
gc(414) = 7
    
gc(442) = 8
gc(443) = 8
gc(441) = 8
gc(444) = 8
    
gc(312) = 9
gc(313) = 9
    
gc(122) = 10
gc(123) = 10
gc(121) = 10
    
gc(221) = 11
gc(224) = 11
gc(322) = 11
gc(323) = 11
gc(321) = 11
gc(324) = 11
    
gc(111) = 12
gc(114) = 12
    
gc(124) = 13
    
gc(222) = 14
gc(223) = 14
    
gc(332) = 15
gc(333) = 15
gc(331) = 15
gc(334) = 15
    
gc(232) = 16
gc(233) = 16
gc(231) = 16
gc(234) = 16
gc(142) = 16
gc(143) = 16
    
gc(132) = 17
gc(133) = 17
gc(131) = 17
gc(134) = 17
    
gc(244) = 18
    
gc(212) = 19
gc(213) = 19
    
gc(422) = 20
gc(423) = 20
gc(421) = 20
gc(424) = 20
   
!
! ------------------------------
! --- aa MW, C and N content ---
! ------------------------------
!
!MWaa(1) = 89 ! not used
!MWaa(2) = 174
!MWaa(3) = 132
!MWaa(4) = 133
!MWaa(5) = 121
!MWaa(6) = 146
!MWaa(7) = 147
!MWaa(8) = 75
!MWaa(9) = 155
!MWaa(10) = 131
!MWaa(11) = 131
!MWaa(12) = 146
!MWaa(13) = 149
!MWaa(14) = 165
!MWaa(15) = 115
!MWaa(16) = 105
!MWaa(17) = 119
!MWaa(18) = 204
!MWaa(19) = 181
!MWaa(20) = 117

raaC(1) = 3.
raaC(2) = 6.
raaC(3) = 4.
raaC(4) = 4.
raaC(5) = 3.
raaC(6) = 5.
raaC(7) = 5.
raaC(8) = 2.
raaC(9) = 6.
raaC(10) = 6.
raaC(11) = 6.
raaC(12) = 6.
raaC(13) = 5.
raaC(14) = 9.
raaC(15) = 5.
raaC(16) = 3.
raaC(17) = 4.
raaC(18) = 11.
raaC(19) = 9.
raaC(20) = 5.

raaN(1) = 1.
raaN(2) = 4.
raaN(3) = 2.
raaN(4) = 1.
raaN(5) = 1.
raaN(6) = 2.
raaN(7) = 1.
raaN(8) = 1.
raaN(9) = 3.
raaN(10) = 1.
raaN(11) = 1.
raaN(12) = 2.
raaN(13) = 1.
raaN(14) = 1.
raaN(15) = 1.
raaN(16) = 1.
raaN(17) = 1.
raaN(18) = 2.
raaN(19) = 1.
raaN(20) = 1.

!reverse aa N req diagnostic simulation
!do iaa1 = 1, 20
!    if (raaN(iaa1).gt.1.) then
!        raaN(iaa1) = 1.
!    else
!        raaN(iaa1) = 2.
!    endif
!enddo

!all same aa N req diagnostic simulation
!do iaa1 = 1, 20
!    raaN(iaa1) = 1.
!enddo

int_saa(1) = 'A'
int_saa(2) = 'R'
int_saa(3) = 'N'
int_saa(4) = 'D'
int_saa(5) = 'C'
int_saa(6) = 'Q'
int_saa(7) = 'E'
int_saa(8) = 'G'
int_saa(9) = 'H'
int_saa(10) = 'I'
int_saa(11) = 'L'
int_saa(12) = 'K'
int_saa(13) = 'M'
int_saa(14) = 'F'
int_saa(15) = 'P'
int_saa(16) = 'S'
int_saa(17) = 'T'
int_saa(18) = 'W'
int_saa(19) = 'Y'
int_saa(20) = 'V'
int_saa(21) = 'START'
int_saa(22) = 'STOP'

int_scod(1) = 'GCT'
int_scod(2) = 'GCC'
int_scod(3) = 'GCA'
int_scod(4) = 'GCG'
int_scod(5) = 'CGT'
int_scod(6) = 'CGC'
int_scod(7) = 'CGA'
int_scod(8) = 'CGG'
int_scod(9) = 'AGA'
int_scod(10) = 'AGG'
int_scod(11) = 'AAT'
int_scod(12) = 'AAC'
int_scod(13) = 'GAT'
int_scod(14) = 'GAC'
int_scod(15) = 'TGT'
int_scod(16) = 'TGC'
int_scod(17) = 'CAA'
int_scod(18) = 'CAG'
int_scod(19) = 'GAA'
int_scod(20) = 'GAG'
int_scod(21) = 'GGT'
int_scod(22) = 'GGC'
int_scod(23) = 'GGA'
int_scod(24) = 'GGG'
int_scod(25) = 'CAT'
int_scod(26) = 'CAC'
int_scod(27) = 'ATT'
int_scod(28) = 'ATC'
int_scod(29) = 'ATA'
int_scod(30) = 'TTA'
int_scod(31) = 'TTG'
int_scod(32) = 'CTT'
int_scod(33) = 'CTC'
int_scod(34) = 'CTA'
int_scod(35) = 'CTG'
int_scod(36) = 'AAA'
int_scod(37) = 'AAG'
int_scod(38) = 'ATG'
int_scod(39) = 'TTT'
int_scod(40) = 'TTC'
int_scod(41) = 'CCT'
int_scod(42) = 'CCC'
int_scod(43) = 'CCA'
int_scod(44) = 'CCG'
int_scod(45) = 'TCT'
int_scod(46) = 'TCC'
int_scod(47) = 'TCA'
int_scod(48) = 'TCG'
int_scod(49) = 'AGT'
int_scod(50) = 'AGC'
int_scod(51) = 'ACT'
int_scod(52) = 'ACC'
int_scod(53) = 'ACA'
int_scod(54) = 'ACG'
int_scod(55) = 'TGG'
int_scod(56) = 'TAT'
int_scod(57) = 'TAC'
int_scod(58) = 'GTT'
int_scod(59) = 'GTC'
int_scod(60) = 'GTA'
int_scod(61) = 'GTG'
int_scod(62) = 'TAA'
int_scod(63) = 'TGA'
int_scod(64) = 'TAG'

int_icod(432) = 1
int_icod(433) = 2
int_icod(431) = 3
int_icod(434) = 4
int_icod(342) = 5
int_icod(343) = 6
int_icod(341) = 7
int_icod(344) = 8
int_icod(141) = 9
int_icod(144) = 10
int_icod(112) = 11
int_icod(113) = 12
int_icod(412) = 13
int_icod(413) = 14
int_icod(242) = 15
int_icod(243) = 16
int_icod(311) = 17
int_icod(314) = 18
int_icod(411) = 19
int_icod(414) = 20
int_icod(442) = 21
int_icod(443) = 22
int_icod(441) = 23
int_icod(444) = 24
int_icod(312) = 25
int_icod(313) = 26
int_icod(122) = 27
int_icod(123) = 28
int_icod(121) = 29
int_icod(221) = 30
int_icod(224) = 31
int_icod(322) = 32
int_icod(323) = 33
int_icod(321) = 34
int_icod(324) = 35
int_icod(111) = 36
int_icod(114) = 37
int_icod(124) = 38
int_icod(222) = 39
int_icod(223) = 40
int_icod(332) = 41
int_icod(333) = 42
int_icod(331) = 43
int_icod(334) = 44
int_icod(232) = 45
int_icod(233) = 46
int_icod(231) = 47
int_icod(234) = 48
int_icod(142) = 49
int_icod(143) = 50
int_icod(132) = 51
int_icod(133) = 52
int_icod(131) = 53
int_icod(134) = 54
int_icod(244) = 55
int_icod(212) = 56
int_icod(213) = 57
int_icod(422) = 58
int_icod(423) = 59
int_icod(421) = 60
int_icod(424) = 61
int_icod(211) = 62
int_icod(241) = 63
int_icod(214) = 64

!
! ------------------------
! --- aa dissimilarity ---
! ------------------------
!
do iaa1 = 1, 22
    do iaa2 = 1, 22
        MdXX(iaa1,iaa2) = 999.
    enddo
enddo

MdXX(1,2) = 2.36
MdXX(1,3) = 1.64
MdXX(1,4) = 2.18
MdXX(1,5) = 1.81
MdXX(1,6) = 1.84
MdXX(1,7) = 2.05
MdXX(1,8) = 1.24
MdXX(1,9) = 2.34
MdXX(1,10) = 2.39
MdXX(1,11) = 1.9
MdXX(1,12) = 2.3
MdXX(1,13) = 1.85
MdXX(1,14) = 2.95
MdXX(1,15) = 0.85
MdXX(1,16) = 0.6
MdXX(1,17) = 1.23
MdXX(1,18) = 4.01
MdXX(1,19) = 2.41
MdXX(1,20) = 1.39
MdXX(2,3) = 2.01
MdXX(2,4) = 2.11
MdXX(2,5) = 2.68
MdXX(2,6) = 1.41
MdXX(2,7) = 1.46
MdXX(2,8) = 2.55
MdXX(2,9) = 1.27
MdXX(2,10) = 2.86
MdXX(2,11) = 2.53
MdXX(2,12) = 1.33
MdXX(2,13) = 2.68
MdXX(2,14) = 3.03
MdXX(2,15) = 2.39
MdXX(2,16) = 2.23
MdXX(2,17) = 2.45
MdXX(2,18) = 3
MdXX(2,19) = 1.84
MdXX(2,20) = 2.47
MdXX(3,4) = 0.8
MdXX(3,5) = 2
MdXX(3,6) = 1.41
MdXX(3,7) = 1.48
MdXX(3,8) = 1.93
MdXX(3,9) = 1.24
MdXX(3,10) = 2.19
MdXX(3,11) = 2.29
MdXX(3,12) = 1.46
MdXX(3,13) = 2.38
MdXX(3,14) = 2.79
MdXX(3,15) = 2.01
MdXX(3,16) = 1.09
MdXX(3,17) = 1.58
MdXX(3,18) = 4.58
MdXX(3,19) = 2.38
MdXX(3,20) = 1.99
MdXX(4,5) = 2.95
MdXX(4,6) = 2.11
MdXX(4,7) = 1.62
MdXX(4,8) = 1.72
MdXX(4,9) = 2.01
MdXX(4,10) = 2.65
MdXX(4,11) = 2.98
MdXX(4,12) = 1.9
MdXX(4,13) = 3.22
MdXX(4,14) = 3.2
MdXX(4,15) = 2.51
MdXX(4,16) = 1.96
MdXX(4,17) = 2.62
MdXX(4,18) = 5.38
MdXX(4,19) = 3.35
MdXX(4,20) = 2.26
MdXX(5,6) = 2.36
MdXX(5,7) = 3.3
MdXX(5,8) = 2.47
MdXX(5,9) = 1.82
MdXX(5,10) = 1.92
MdXX(5,11) = 1.95
MdXX(5,12) = 2.87
MdXX(5,13) = 2.1
MdXX(5,14) = 2.23
MdXX(5,15) = 1.74
MdXX(5,16) = 1.21
MdXX(5,17) = 1.67
MdXX(5,18) = 3.99
MdXX(5,19) = 1.7
MdXX(5,20) = 1.78
MdXX(6,7) = 1.04
MdXX(6,8) = 2.73
MdXX(6,9) = 1.25
MdXX(6,10) = 2.75
MdXX(6,11) = 2.16
MdXX(6,12) = 1.31
MdXX(6,13) = 2.05
MdXX(6,14) = 3.22
MdXX(6,15) = 2.23
MdXX(6,16) = 1.54
MdXX(6,17) = 1.52
MdXX(6,18) = 2.94
MdXX(6,19) = 1.99
MdXX(6,20) = 2.32
MdXX(7,8) = 2.41
MdXX(7,9) = 2
MdXX(7,10) = 3.25
MdXX(7,11) = 2.82
MdXX(7,12) = 1.12
MdXX(7,13) = 2.65
MdXX(7,14) = 3.73
MdXX(7,15) = 2.64
MdXX(7,16) = 2.18
MdXX(7,17) = 2.59
MdXX(7,18) = 3.79
MdXX(7,19) = 2.79
MdXX(7,20) = 2.58
MdXX(8,9) = 2.74
MdXX(8,10) = 3.06
MdXX(8,11) = 3.02
MdXX(8,12) = 3.1
MdXX(8,13) = 3.2
MdXX(8,14) = 3.3
MdXX(8,15) = 1.77
MdXX(8,16) = 1.4
MdXX(8,17) = 2.57
MdXX(8,18) = 5.05
MdXX(8,19) = 3.43
MdXX(8,20) = 2.06
MdXX(9,10) = 2.49
MdXX(9,11) = 2.31
MdXX(9,12) = 1.73
MdXX(9,13) = 2.55
MdXX(9,14) = 2.62
MdXX(9,15) = 2.58
MdXX(9,16) = 1.77
MdXX(9,17) = 1.8
MdXX(9,18) = 3.49
MdXX(9,19) = 1.36
MdXX(9,20) = 2.4
MdXX(10,11) = 0.47
MdXX(10,12) = 2.24
MdXX(10,13) = 1.12
MdXX(10,14) = 0.71
MdXX(10,15) = 2.9
MdXX(10,16) = 2.41
MdXX(10,17) = 2.45
MdXX(10,18) = 3.9
MdXX(10,19) = 1.53
MdXX(10,20) = 0.82
MdXX(11,12) = 2
MdXX(11,13) = 0.47
MdXX(11,14) = 1.23
MdXX(11,15) = 2.69
MdXX(11,16) = 2.1
MdXX(11,17) = 1.79
MdXX(11,18) = 3.05
MdXX(11,19) = 1.17
MdXX(11,20) = 0.75
MdXX(12,13) = 1.77
MdXX(12,14) = 2.9
MdXX(12,15) = 2.68
MdXX(12,16) = 2.18
MdXX(12,17) = 2.09
MdXX(12,18) = 3.89
MdXX(12,19) = 2.26
MdXX(12,20) = 1.99
MdXX(13,14) = 1.89
MdXX(13,15) = 2.48
MdXX(13,16) = 2.04
MdXX(13,17) = 1.46
MdXX(13,18) = 3.05
MdXX(13,19) = 1.43
MdXX(13,20) = 0.94
MdXX(14,15) = 3.6
MdXX(14,16) = 2.78
MdXX(14,17) = 2.93
MdXX(14,18) = 3.69
MdXX(14,19) = 1.47
MdXX(14,20) = 1.47
MdXX(15,16) = 0.99
MdXX(15,17) = 1.46
MdXX(15,18) = 4.33
MdXX(15,19) = 2.88
MdXX(15,20) = 2.21
MdXX(16,17) = 0.82
MdXX(16,18) = 4.07
MdXX(16,19) = 2.06
MdXX(16,20) = 1.88
MdXX(17,18) = 3.6
MdXX(17,19) = 2.06
MdXX(17,20) = 1.87
MdXX(18,19) = 2.04
MdXX(18,20) = 3.92
MdXX(19,20) = 1.95

do iaa1 = 1, 19
    do iaa2 = iaa1 + 1, 20
        MdXX(iaa2,iaa1) = MdXX(iaa1,iaa2)
    enddo
enddo

!
! ---------------------
! --- serial number ---
! ---------------------
!
gsn = 1
gsnmax = 999999999

!
! ----------------
! --- mutation ---
! ----------------
!
do nt1 = 1, 4
    do nt2 = 1, 4
        uEbij(nt1,nt2) = 0.
        if (nt1.ne.nt2) then
            uEbij(nt1,nt2) = ubij(nt1,nt2) * 2. * fu 
        endif
    enddo
enddo
uEbk(1) = uEbij(1,2) + uEbij(1,3) + uEbij(1,4)
uEbk(2) = uEbij(2,1) + uEbij(2,3) + uEbij(2,4)
uEbk(3) = uEbij(3,1) + uEbij(3,2) + uEbij(3,4)
uEbk(4) = uEbij(4,1) + uEbij(4,2) + uEbij(4,3)

!
! ---------------------
! --- recombination ---
! ---------------------
!
uREg = uRg * 2.

!
! --- changes ---
!
ncg_sp_mx = 0
icg = 0

!
! ------------------------------
! --- DNA C, N and P content ---
! ------------------------------
!
NA = 6.022e23

rGCC = 19. * n0DNA / NA
rATC = 20. * n0DNA / NA
rGCN = 8. * n0DNA / NA
rATN = 7. * n0DNA / NA
rGCP = 2. * n0DNA / NA
rATP = 2. * n0DNA / NA

!
! ----------------
! --- dilution ---
! ----------------
!
tDNext = dtD
ipD = 1



!
! --- extracellular concentrations ---
!
do isp1 = 1, dimnsp
    SC_sp(isp1) = SC0
    SN_sp(isp1) = SN0
    SP_sp(isp1) = SP0
enddo

!
! --- code performance ---
!
cp_delta = 0.
cp_init_delta = 0.
cp_iter_delta = 0.
cp_fix_delta = 0.
cp_output_delta = 0.
cp_genomes_delta = 0.
cp_partot_delta = 0.
cp_parwait_delta = 0.
cp_mix1_delta = 0.
cp_mix2_delta = 0.
do isp1 = 1, dimnsp
    cp_ENS_delta_sp(isp1) = 0.
    cp_HANS1_delta_sp(isp1) = 0.
    cp_HANS2_delta_sp(isp1) = 0.
    !cp_HANS2a_delta_sp(isp1) = 0.
    !cp_HANS2b_delta_sp(isp1) = 0.
    !cp_HANS2c_delta_sp(isp1) = 0.
    cp_mr1_delta_sp(isp1) = 0.
    cp_mr2_delta_sp(isp1) = 0.
    !cp_mr2a_delta_sp(isp1) = 0.
    !cp_mr2b_delta_sp(isp1) = 0.
    !cp_mr2c_delta_sp(isp1) = 0.
enddo

cp_t1 = omp_get_wtime()
cp_init_t1 = omp_get_wtime()

!
! --- sarkar preset init ---
!
if (oHANSs.eq.3) then
    isarkar = 0
else
    isarkar = 1
endif

!
! ---------------------------
! ---------------------------
! --- time-variable input ---
! ---------------------------
! ---------------------------
!
open(112,file='IAM_I_t.txt')
it = 1
do while (.not.eof(112))
    if (it.gt.dimntb) then
        write(*,*) 'IAM\T: ERR: it > dimntb.', it, dimntb
        goto 9100
    endif
    read(112,*) tt(it), &
        tfQ(it), tfSinC(it), tfSinN(it), tfSinP(it), &
        tfSCfix(it), tfSNfix(it), tfSPfix(it), &
        tfX1(it), tfX2(it)
    !write(*,*) 'IAM\T: it, tt = ', it, tt(it)
    it = it + 1
enddo
close(112)

it = 1

!
! ------------------
! ------------------
! --- no. agents ---
! ------------------
! ------------------
!
nic = Popic
if ((oHANS.eq.1).and.(oStartDNA.gt.1)) then
    nic = 10 * dimnsp
endif

!
! --------------------------------------
! --------------------------------------
! --- initial proteins and genome(s) ---
! --------------------------------------
! --------------------------------------
!
write(*,*) 'IAM\INIT:'
write(*,*) 'IAM\INIT: Initial proteins and genome(s).'

ngt = 0
if (oStartDNA.le.3) then
    write(*,*) 'IAM\INIT:'
    write(*,*) 'IAM\INIT: Generating synthetic.'
    nnt = dimnnt
    if (oStartDNA.eq.1) then
        ngt = nic
    else
        ngt = 1
    endif
    do ig1 = 1, ngt
        if ((oStartDNA.eq.1).or.(oStartDNA.eq.2)) then
            do int1 = 1, nnt
                r = r4_uni(idum)
                if (r.lt.ficA) then
                    nt1 = 1
                elseif (r.lt.(ficA+ficT)) then
                    nt1 = 2
                elseif (r.lt.(ficA+ficT+ficC)) then
                    nt1 = 3
                else
                    nt1 = 4
                endif
                g_nt(ig1,int1) = nt1
            enddo
        else
            nt1 = 1
            do int1 = 1, nnt
                g_nt(ig1,int1) = nt1
            enddo
        endif
        ip1 = 1
        int1 = 1
        do while(.true.)
            r = r4_uni(idum)
            if (r.lt..5) then
                g_pother(ig1,ip1) = 1
            else
                g_pother(ig1,ip1) = -1
            endif
            r = r4_uni(idum)
            g_pgamma(ig1,ip1) = 2. * r
            r = r4_uni(idum)
            pspacing = pminspacing + (pmaxspacing - pminspacing) * r
            g_pntstart(ig1,ip1) = int1 + pspacing + 1
            if (g_pntstart(ig1,ip1).gt.nnt) then
                g_np(ig1) = ip1 - 1
                goto 4030
            endif
            r = r4_uni(idum)
            plength = nint((pminlength + (pmaxlength - pminlength) * r) / 3.) * 3
            g_pntstop(ig1,ip1) = g_pntstart(ig1,ip1) + plength - 1
            if (g_pntstop(ig1,ip1).gt.nnt) then
                g_np(ig1) = ip1 - 1
                goto 4030
            endif
            int1 = g_pntstop(ig1,ip1)
            ip1 = ip1 + 1
        enddo
4030 continue
        if (g_np(ig1).eq.0) then
            write(*,*) 'IAM\INIT: ERR. g_np = 0. ig1 = ', ig1
            goto 9100
        endif
    enddo ! ig1
else
    write(*,*) 'IAM\INIT:'
    write(*,*) 'IAM\INIT: Reading from file.'
    ngt = ngt + 1
    ig1 = ngt
    nnt = 0
    ip2 = 0
    do iger = 1, nger
        write(*,*) 'IAM\INIT:'
        write(*,*) 'IAM\INIT: iger = ', iger
        open(122,file='IAM_I_r_'//sap(iger)//'.txt')
        do while(not(eof(122)))
            read(122,*) ip1, pntstartr, pntstopr, potherr, pgammar
            !write(*,*) ip1, pntstartr, pntstopr, potherr, pgammar
            ip2 = ip2 + 1
            int1 = nnt + pntstartr
            int2 = nnt + pntstopr
            if (int1.ge.dimnnt) then
                ip2 = ip2 - 1
                write(*,*) 'IAM\INIT: Warning: Extra protein(s)'
                write(*,*) 'IAM\INIT: int1, int2, dimnnt = ', int1, int2, dimnnt
                goto 4020
            endif
            if (int2.gt.dimnnt) then
                write(*,*) 'IAM\INIT: Warning: Truncating protein'
                write(*,*) 'IAM\INIT: int1, int2, dimnnt = ', int1, int2, dimnnt
                int2 = dimnnt
            endif
            g_pntstart(ig1,ip2) = int1
            g_pntstop(ig1,ip2) = int2
            g_pother(ig1,ip2) = potherr
            int3 = g_pntstop(ig1,ip2) - g_pntstart(ig1,ip2) + 1
            if (((dble(int3)/3.)-dble(int3/3)).ne.0.) then
                write(*,*) 'IAM\INIT: ERR: Protein not multiple of 3.'
                write(*,*) 'IAM\INIT: ', ip2, g_pntstart(ig1,ip2), g_pntstop(ig1,ip2)
                write(*,*) 'IAM\INIT: ', dble(int3)/3., dble(int3/3)
                goto 9100
            endif
            g_pgamma(ig1,ip2) = pgammar
            if (ip1.eq.1) then
                write(*,*) 'IAM\INIT: First protein start, stop, other, gamma:'
                write(*,*) 'IAM\INIT: read  = ', pntstartr, pntstopr, potherr, pgammar
                write(*,*) 'IAM\INIT: corr  = ', g_pntstart(ig1,ip2), g_pntstop(ig1,ip2), g_pother(ig1,ip2), g_pgamma(ig1,ip2)
            endif
        enddo
4020 continue
        close(122)
        write(*,*) 'IAM\INIT: last protein start, stop, other, gamma:'
        write(*,*) 'IAM\INIT: read  = ', pntstartr, pntstopr, potherr, pgammar
        write(*,*) 'IAM\INIT: corr  = ', g_pntstart(ig1,ip2), g_pntstop(ig1,ip2), g_pother(ig1,ip2), g_pgamma(ig1,ip2)
        open(121,file='IAM_I_g_'//sap(iger)//'.txt')
        read(121,208) shr
        write(*,*) 'IAM\INIT: Genome header:'
        write(*,*) shr
        do while(.true.)
            read(121,207) (sntr(int1),int1=1,60)
            !write(*,*) (sntr(int1),int1=1,60)
            do int1 = 1, 60
                if (sntr(int1).eq.'') then
                    write(*,*) 'IAM\INIT: Stop on blank.'
                    int2 = int1 - 1
                    goto 4050
                endif
                nnt = nnt + 1
                if (sntr(int1).eq.'A') then
                    g_nt(ig1,nnt) = 1
                elseif (sntr(int1).eq.'T') then
                    g_nt(ig1,nnt) = 2
                elseif (sntr(int1).eq.'C') then
                    g_nt(ig1,nnt) = 3
                elseif (sntr(int1).eq.'G') then
                    g_nt(ig1,nnt) = 4
                else
                    write(*,*) 'IAM\INIT: ERR. Non-ATCG in genome: ', nnt, sntr(int1)
                    goto 9100
                endif
                if (nnt.eq.dimnnt) then
                    write(*,*) 'IAM\INIT: Stop on dimnnt.'
                    int2 = int1
                    goto 4050
                endif
            enddo ! int1
            if (eof(121)) then
                write(*,*) 'IAM\INIT: Stop on eof.'
                int2 = int1
                goto 4050
            endif
        enddo
4050 continue
        close(121)
        write(*,*) 'IAM\INIT: Last few nt read  = ', (sntr(int3),int3=max(1,int2-9),int2)
        write(*,*) 'IAM\INIT: Last 10 nt in mem = ', (int_snt(g_nt(ig1,int3)),int3=nnt-9,nnt)
        write(*,*) 'IAM\INIT: total bp read = ', nnt
        if (nnt.eq.dimnnt) then
            goto 4060
        endif
    enddo ! iger
4060 continue
    g_np(ig1) = ip2
    write(*,*) 'IAM\INIT: bp processed = ', nnt
    if (nnt.ne.dimnnt) then
        write(*,*) "IAM\INIT: ERR: Does not match dimnnt = ", dimnnt
        goto 9100
    endif
endif

write(*,*) 'IAM\INIT: Initial genomes done.'
write(*,*) 'IAM\INIT: ngt = ', ngt

!
! --- perturbation analysis ---
!
! notea:
! - only set up to work with initial population and one genome
! - use jdum/RSj, which is same for all runs, so same changes are made
!
if (pfAT.gt.0.) then
    write(*,*) 'IAM\PER:'
    write(*,*) 'IAM\PER: pfAT = ', pfAT
    ig1 = 1
    do int1 = 1, dimnnt
        nt1 = g_nt(ig1,int1)
        if (nt1.le.2) then
            r = r4_uni(jdum)
            if (r.lt.pfAT) then
                r = r4_uni(jdum)
                if (r.lt.0.5) then
                    nt2 = 3
                else
                    nt2 = 4
                endif
                g_nt(ig1,int1) = nt2
            endif
        endif
    enddo
endif    
if (pfGC.gt.0.) then
    write(*,*) 'IAM\PER:'
    write(*,*) 'IAM\PER: pfGC = ', pfGC
    ig1 = 1
    do int1 = 1, dimnnt
        nt1 = g_nt(ig1,int1)
        if (nt1.ge.3) then
            r = r4_uni(jdum)
            if (r.lt.pfGC) then
                r = r4_uni(jdum)
                if (r.lt.0.5) then
                    nt2 = 1
                else
                    nt2 = 2
                endif
                g_nt(ig1,int1) = nt2
            endif
        endif
    enddo
endif    

!
! --- process genome(s) ---
!
write(*,*) 'IAM\GEN:'
write(*,*) 'IAM\GEN: *** Genome properties ***'
write(*,*) 'IAM\GEN:'
uEgmin = 1.e9
uEgmax = -1.e9
do ig1 = 1, ngt
    g_ntn(ig1) = 0.
    g_ntnS(ig1) = 0.
    g_ntnN(ig1) = 0.
    g_ntnX(ig1) = 0.
    do nt1 = 1, 4
        g_ntnk(ig1,nt1) = 0.
        g_ntnSk(ig1,nt1) = 0.
        g_ntnNk(ig1,nt1) = 0.
        g_ntnXk(ig1,nt1) = 0.
    enddo
    do nt1 = 1, 4
        do nt2 = 1, 4
            g_ntnij(ig1,nt1,nt2) = 0.
            g_ntnSij(ig1,nt1,nt2) = 0.
            g_ntnNij(ig1,nt1,nt2) = 0.
            g_ntnXij(ig1,nt1,nt2) = 0.
        enddo
    enddo
    g_ncf(ig1) = 0
    g_ncff(ig1) = 0
    g_ncSff(ig1) = 0
    g_ncNff(ig1) = 0
    g_ncmff(ig1) = 0
    g_ncrff(ig1) = 0
    g_uEg(ig1) = 0.
    g_q0DNAC(ig1) = 0.
    g_q0DNAN(ig1) = 0.
    g_q0DNAP(ig1) = 0.
    g_pntn(ig1,:) = 0.
    g_pntnS(ig1,:) = 0.
    g_pntnN(ig1,:) = 0.
    do int1 = 1, dimnnt
        g_nt_o(ig1,int1) = g_nt(ig1,int1)
        g_ntfS(ig1,int1) = 0
        g_cmr(ig1,int1) = 0
        nt1 = g_nt(ig1,int1)
        g_uEg(ig1) = g_uEg(ig1) + uEbk(nt1)
        ! no. total, S, N and X sites
        g_ntn(ig1) = g_ntn(ig1) + 1.
        g_ntnk(ig1,nt1) = g_ntnk(ig1,nt1) + 1.
        do int3 = 1, 3
            nt2 = cint(nt1,int3)
            g_ntnij(ig1,nt1,nt2) = g_ntnij(ig1,nt1,nt2) + 1.
        enddo
        g_ip(ig1,int1) = -9
        do ip1 = 1, g_np(ig1)
            if ((int1.ge.g_pntstart(ig1,ip1)).and.(int1.le.g_pntstop(ig1,ip1))) then ! protein found
                g_ip(ig1,int1) = ip1
                if (g_pother(ig1,ip1).eq.-1) then
                    int2 = g_pntstart(ig1,ip1)+(int1-g_pntstart(ig1,ip1))/3*3.
                    if (int1.eq.int2) then
                        nta(1) = nt1
                        nta(2) = g_nt(ig1,int1+1)
                        nta(3) = g_nt(ig1,int1+2)
                        int4 = 1
                    elseif (int1.eq.int2+1) then
                        nta(1) = g_nt(ig1,int1-1)
                        nta(2) = nt1
                        nta(3) = g_nt(ig1,int1+1)
                        int4 = 2
                    else
                        nta(1) = g_nt(ig1,int1-2)
                        nta(2) = g_nt(ig1,int1-1)
                        nta(3) = nt1
                        int4 = 3
                    endif
                    cod1 = nta(1)*100+nta(2)*10+nta(3)*1
                else
                    int2 = g_pntstop(ig1,ip1)-(g_pntstop(ig1,ip1)-int1)/3*3.
                    if (int1.eq.int2) then
                        nta(1) = nt1
                        nta(2) = g_nt(ig1,int1-1)
                        nta(3) = g_nt(ig1,int1-2)
                        int4 = 1
                    elseif (int1.eq.int2-1) then
                        nta(1) = g_nt(ig1,int1+1)
                        nta(2) = nt1
                        nta(3) = g_nt(ig1,int1-1)
                        int4 = 2
                    else
                        nta(1) = g_nt(ig1,int1+2)
                        nta(2) = g_nt(ig1,int1+1)
                        nta(3) = nt1
                        int4 = 3
                    endif
                    cod1 = int_opo(nta(1))*100+int_opo(nta(2))*10+int_opo(nta(3))*1
                endif
                iaa1 = gc(cod1)
                do int3 = 1, 3
                    nt2 = cint(nt1,int3)
                    if (int4.eq.1) then
                        ntb(1) = nt2
                        ntb(2) = nta(2)
                        ntb(3) = nta(3)
                    elseif (int4.eq.2) then
                        ntb(1) = nta(1)
                        ntb(2) = nt2
                        ntb(3) = nta(3)
                    else
                        ntb(1) = nta(1)
                        ntb(2) = nta(2)
                        ntb(3) = nt2
                    endif
                    if (g_pother(ig1,ip1).eq.-1) then
                        cod2 = ntb(1)*100+ntb(2)*10+ntb(3)*1
                    else
                        cod2 = int_opo(ntb(1))*100+int_opo(ntb(2))*10+int_opo(ntb(3))*1
                    endif
                    iaa2 = gc(cod2)
                    if (iaa1.eq.iaa2) then
                        !write(*,*) 'IAM\GEN: Synonymous.'
                        g_ntfS(ig1,int1) = g_ntfS(ig1,int1) + 1
                        g_ntnS(ig1) = g_ntnS(ig1) + 1./3.
                        g_ntnSk(ig1,nt1) = g_ntnSk(ig1,nt1) + 1./3.
                        g_ntnSij(ig1,nt1,nt2) = g_ntnSij(ig1,nt1,nt2) + 1.
                        g_pntnS(ig1,ip1) = g_pntnS(ig1,ip1) + 1./3.
                    else
                        !write(*,*) 'IAM\GEN: Nonsynonymous.'
                        g_ntnN(ig1) = g_ntnN(ig1) + 1./3.
                        g_ntnNk(ig1,nt1) = g_ntnNk(ig1,nt1) + 1./3.
                        g_ntnNij(ig1,nt1,nt2) = g_ntnNij(ig1,nt1,nt2) + 1.
                        g_pntnN(ig1,ip1) = g_pntnN(ig1,ip1) + 1./3.
                    endif
                enddo ! int3
                goto 4075
            endif ! protein found
        enddo ! ip1
4075 continue
        if (nt1.le.2) then
            g_q0DNAC(ig1) = g_q0DNAC(ig1) + rATC
            g_q0DNAN(ig1) = g_q0DNAN(ig1) + rATN
            g_q0DNAP(ig1) = g_q0DNAP(ig1) + rATP
        else
            g_q0DNAC(ig1) = g_q0DNAC(ig1) + rGCC
            g_q0DNAN(ig1) = g_q0DNAN(ig1) + rGCN
            g_q0DNAP(ig1) = g_q0DNAP(ig1) + rGCP
        endif
    enddo ! int1
    g_ntnX(ig1) = g_ntn(ig1) - g_ntnS(ig1) - g_ntnN(ig1)
    do nt1 = 1, 4
        g_ntnXk(ig1,nt1) = g_ntnk(ig1,nt1) - g_ntnSk(ig1,nt1) - g_ntnNk(ig1,nt1)
    enddo
    do nt1 = 1, 4
        do nt2 = 1, 4
            if (nt1.ne.nt2) then
                g_ntnXij(ig1,nt1,nt2) = g_ntnij(ig1,nt1,nt2) - g_ntnSij(ig1,nt1,nt2) - g_ntnNij(ig1,nt1,nt2)
            endif
        enddo
    enddo
    if (g_uEg(ig1).lt.uEgmin) uEgmin = g_uEg(ig1)
    if (g_uEg(ig1).gt.uEgmax) uEgmax = g_uEg(ig1)
    write(*,*) 'IAM\GEN:'
    write(*,*) 'IAM\GEN: ig1 = ', ig1
    write(*,*) 'IAM\GEN:'
    write(*,*) 'IAM\GEN: * mutation and recombination *'
    write(*,*) 'IAM\GEN:'
    write(*,*) 'IAM\GEN: uEg =        ', g_uEg(ig1)
    write(*,*) 'IAM\GEN: uEg*fM =     ', g_uEg(ig1)*fM
    write(*,*) 'IAM\GEN: uEg*fM*xM =  ', g_uEg(ig1)*fM*xM
    write(*,*) 'IAM\GEN: uREg =       ', uREg
    write(*,*) 'IAM\GEN: uREg*fR =    ', uREg*fR
    write(*,*) 'IAM\GEN: uREg*fR*xR = ', uREg*fR*xR
    write(*,*) 'IAM\GEN:'
    write(*,*) 'IAM\GEN: * ATCG composition *'
    write(*,*) 'IAM\GEN:'
    write(*,*) 'IAM\GEN: Genome:'
    write(*,*) 'IAM\GEN: n, f A =   ', g_ntnk(ig1,1), dble(g_ntnk(ig1,1)) / dble(dimnnt)
    write(*,*) 'IAM\GEN: n, f T =   ', g_ntnk(ig1,2), dble(g_ntnk(ig1,2)) / dble(dimnnt)
    write(*,*) 'IAM\GEN: n, f C =   ', g_ntnk(ig1,3), dble(g_ntnk(ig1,3)) / dble(dimnnt)
    write(*,*) 'IAM\GEN: n, f G =   ', g_ntnk(ig1,4), dble(g_ntnk(ig1,4)) / dble(dimnnt)
    write(*,*) 'IAM\GEN: n, f GC =  ', (g_ntnk(ig1,3)+g_ntnk(ig1,4)), dble((g_ntnk(ig1,3)+g_ntnk(ig1,4))) / dble(dimnnt)
    write(*,*) 'IAM\GEN:'
    write(*,*) 'IAM\GEN: by TXSN:'
    write(*,*) 'IAM\GEN: T'
    write(*,*) 'IAM\GEN: n, f A =   ', g_ntnk(ig1,1), g_ntnk(ig1,1) / g_ntn(ig1)
    write(*,*) 'IAM\GEN: n, f T =   ', g_ntnk(ig1,2), g_ntnk(ig1,2) / g_ntn(ig1)
    write(*,*) 'IAM\GEN: n, f C =   ', g_ntnk(ig1,3), g_ntnk(ig1,3) / g_ntn(ig1)
    write(*,*) 'IAM\GEN: n, f G =   ', g_ntnk(ig1,4), g_ntnk(ig1,4) / g_ntn(ig1)
    write(*,*) 'IAM\GEN: n, f GC =  ', (g_ntnk(ig1,3)+g_ntnk(ig1,4)), (g_ntnk(ig1,3)+g_ntnk(ig1,4)) / g_ntn(ig1)
    write(*,*) 'IAM\GEN: X'
    write(*,*) 'IAM\GEN: n, f A =   ', g_ntnXk(ig1,1), g_ntnXk(ig1,1) / g_ntnX(ig1)
    write(*,*) 'IAM\GEN: n, f T =   ', g_ntnXk(ig1,2), g_ntnXk(ig1,2) / g_ntnX(ig1)
    write(*,*) 'IAM\GEN: n, f C =   ', g_ntnXk(ig1,3), g_ntnXk(ig1,3) / g_ntnX(ig1)
    write(*,*) 'IAM\GEN: n, f G =   ', g_ntnXk(ig1,4), g_ntnXk(ig1,4) / g_ntnX(ig1)
    write(*,*) 'IAM\GEN: n, f GC =  ', (g_ntnXk(ig1,3)+g_ntnXk(ig1,4)), (g_ntnXk(ig1,3)+g_ntnXk(ig1,4)) / g_ntnX(ig1)
    write(*,*) 'IAM\GEN: S'
    write(*,*) 'IAM\GEN: n, f A =   ', g_ntnSk(ig1,1), g_ntnSk(ig1,1) / g_ntnS(ig1)
    write(*,*) 'IAM\GEN: n, f T =   ', g_ntnSk(ig1,2), g_ntnSk(ig1,2) / g_ntnS(ig1)
    write(*,*) 'IAM\GEN: n, f C =   ', g_ntnSk(ig1,3), g_ntnSk(ig1,3) / g_ntnS(ig1)
    write(*,*) 'IAM\GEN: n, f G =   ', g_ntnSk(ig1,4), g_ntnSk(ig1,4) / g_ntnS(ig1)
    write(*,*) 'IAM\GEN: n, f GC =  ', (g_ntnSk(ig1,3)+g_ntnSk(ig1,4)), (g_ntnSk(ig1,3)+g_ntnSk(ig1,4)) / g_ntnS(ig1)
    write(*,*) 'IAM\GEN: XS'
    write(*,*) 'IAM\GEN: n, f A =   ', g_ntnXk(ig1,1)+g_ntnSk(ig1,1), (g_ntnXk(ig1,1)+g_ntnSk(ig1,1)) / (g_ntnX(ig1)+g_ntnS(ig1))
    write(*,*) 'IAM\GEN: n, f T =   ', g_ntnXk(ig1,2)+g_ntnSk(ig1,2), (g_ntnXk(ig1,2)+g_ntnSk(ig1,2)) / (g_ntnX(ig1)+g_ntnS(ig1))
    write(*,*) 'IAM\GEN: n, f C =   ', g_ntnXk(ig1,3)+g_ntnSk(ig1,3), (g_ntnXk(ig1,3)+g_ntnSk(ig1,3)) / (g_ntnX(ig1)+g_ntnS(ig1))
    write(*,*) 'IAM\GEN: n, f G =   ', g_ntnXk(ig1,4)+g_ntnSk(ig1,4), (g_ntnXk(ig1,4)+g_ntnSk(ig1,4)) / (g_ntnX(ig1)+g_ntnS(ig1))
    write(*,*) 'IAM\GEN: n, f GC =  ', (g_ntnXk(ig1,3)+g_ntnXk(ig1,4))+(g_ntnSk(ig1,3)+g_ntnSk(ig1,4)), ((g_ntnXk(ig1,3)+g_ntnXk(ig1,4))+(g_ntnSk(ig1,3)+g_ntnSk(ig1,4))) / (g_ntnX(ig1)+g_ntnS(ig1))
    write(*,*) 'IAM\GEN: N'
    write(*,*) 'IAM\GEN: n, f A =   ', g_ntnNk(ig1,1), g_ntnNk(ig1,1) / g_ntnN(ig1)
    write(*,*) 'IAM\GEN: n, f T =   ', g_ntnNk(ig1,2), g_ntnNk(ig1,2) / g_ntnN(ig1)
    write(*,*) 'IAM\GEN: n, f C =   ', g_ntnNk(ig1,3), g_ntnNk(ig1,3) / g_ntnN(ig1)
    write(*,*) 'IAM\GEN: n, f G =   ', g_ntnNk(ig1,4), g_ntnNk(ig1,4) / g_ntnN(ig1)
    write(*,*) 'IAM\GEN: n, f GC =  ', (g_ntnNk(ig1,3)+g_ntnNk(ig1,4)), (g_ntnNk(ig1,3)+g_ntnNk(ig1,4)) / g_ntnN(ig1)
    write(*,*) 'IAM\GEN:'
    write(*,*) 'IAM\GEN: * TXSN composition *'
    write(*,*) 'IAM\GEN:'
    write(*,*) 'IAM\GEN: Genome:'
    write(*,*) 'IAM\GEN: n, f T =   ', g_ntn(ig1), g_ntn(ig1)/g_ntn(ig1)
    write(*,*) 'IAM\GEN: n, f X =   ', g_ntnX(ig1), g_ntnX(ig1)/g_ntn(ig1)
    write(*,*) 'IAM\GEN: n, f S =   ', g_ntnS(ig1), g_ntnS(ig1)/g_ntn(ig1)
    write(*,*) 'IAM\GEN: n, f XS =  ', g_ntnX(ig1)+g_ntnS(ig1), (g_ntnX(ig1)+g_ntnS(ig1))/g_ntn(ig1)
    write(*,*) 'IAM\GEN: n, f N =   ', g_ntnN(ig1), g_ntnN(ig1)/g_ntn(ig1)
    write(*,*) 'IAM\GEN:'
    write(*,*) 'IAM\GEN: by ATCG:'
    do nt1 = 1, 4
        write(*,*) 'IAM\GEN: nt1 = ', nt1
        write(*,*) 'IAM\GEN: n, f T =   ', g_ntnk(ig1,nt1), g_ntnk(ig1,nt1)/g_ntnk(ig1,nt1)
        write(*,*) 'IAM\GEN: n, f X =   ', g_ntnXk(ig1,nt1), g_ntnXk(ig1,nt1)/g_ntnk(ig1,nt1)
        write(*,*) 'IAM\GEN: n, f S =   ', g_ntnSk(ig1,nt1), g_ntnSk(ig1,nt1)/g_ntnk(ig1,nt1)
        write(*,*) 'IAM\GEN: n, f XS =  ', g_ntnXk(ig1,nt1)+g_ntnSk(ig1,nt1), (g_ntnXk(ig1,nt1)+g_ntnSk(ig1,nt1))/g_ntnk(ig1,nt1)
        write(*,*) 'IAM\GEN: n, f N =   ', g_ntnNk(ig1,nt1), g_ntnNk(ig1,nt1)/g_ntnk(ig1,nt1)
    enddo
!    write(*,*) 'IAM\GEN: AT fS = ', (g_ntnSk(ig1,1)+g_ntnSk(ig1,2))/(g_ntnk(ig1,1)+g_ntnk(ig1,2))
!    write(*,*) 'IAM\GEN: GC fS = ', (g_ntnSk(ig1,3)+g_ntnSk(ig1,4))/(g_ntnk(ig1,3)+g_ntnk(ig1,4))
    do nt1 = 1, 4
        do nt2 = 1, 4
            if (nt1.ne.nt2) then
                write(*,*) 'IAM\GEN: nt1, nt2 = ', nt1, nt2
                write(*,*) 'IAM\GEN: T n, f =   ', g_ntnij(ig1,nt1,nt2), g_ntnij(ig1,nt1,nt2)/g_ntnij(ig1,nt1,nt2)
                write(*,*) 'IAM\GEN: X n, f =   ', g_ntnXij(ig1,nt1,nt2), g_ntnXij(ig1,nt1,nt2)/g_ntnij(ig1,nt1,nt2)
                write(*,*) 'IAM\GEN: S n, f =   ', g_ntnSij(ig1,nt1,nt2), g_ntnSij(ig1,nt1,nt2)/g_ntnij(ig1,nt1,nt2)
                write(*,*) 'IAM\GEN: XS n, f =  ', g_ntnXij(ig1,nt1,nt2)+g_ntnSij(ig1,nt1,nt2), (g_ntnXij(ig1,nt1,nt2)+g_ntnSij(ig1,nt1,nt2))/g_ntnij(ig1,nt1,nt2)
                write(*,*) 'IAM\GEN: N n, f =   ', g_ntnNij(ig1,nt1,nt2), g_ntnNij(ig1,nt1,nt2)/g_ntnij(ig1,nt1,nt2)
            endif
        enddo
    enddo
    write(*,*) 'IAM\GEN:'
    write(*,*) 'IAM\GEN: * Site x mutation bias *'
    write(*,*) 'IAM\GEN: AT>GC'
    int1 = 1
    nt11 = 1
    nt21 = 4
    nt12 = 2
    nt22 = 3
    bS(int1) = (g_ntnSij(ig1,nt11,nt21)+g_ntnSij(ig1,nt12,nt22))
    bN(int1) = (g_ntnNij(ig1,nt11,nt21)+g_ntnNij(ig1,nt12,nt22))
    bX(int1) = (g_ntnXij(ig1,nt11,nt21)+g_ntnXij(ig1,nt12,nt22))
    bSX(int1) = bS(int1) + bX(int1)
    bT(int1) = bN(int1) + bSX(int1)
    bSxu(int1) = bS(int1) * ubij(nt11,nt21)
    bNxu(int1) = bN(int1) * ubij(nt11,nt21)
    bXxu(int1) = bX(int1) * ubij(nt11,nt21)
    bSXxu(int1) = bSxu(int1) + bXxu(int1)
    bTxu(int1) = bNxu(int1) + bSXxu(int1)
    write(*,*) 'IAM\GEN: bS, bS x ubij   = ', bS(int1), bSxu(int1)
    write(*,*) 'IAM\GEN: bN, bN x ubij   = ', bN(int1), bNxu(int1)
    write(*,*) 'IAM\GEN: bX, bX x ubij   = ', bX(int1), bXxu(int1)
    write(*,*) 'IAM\GEN: bSX, bSX x ubij = ', bSX(int1), bSXxu(int1)
    write(*,*) 'IAM\GEN: bT, bT x ubij   = ', bT(int1), bTxu(int1)
    write(*,*) 'IAM\GEN: GC>AT'
    int1 = 2
    nt11 = 4
    nt21 = 1
    nt12 = 3
    nt22 = 2
    bS(int1) = (g_ntnSij(ig1,nt11,nt21)+g_ntnSij(ig1,nt12,nt22))
    bN(int1) = (g_ntnNij(ig1,nt11,nt21)+g_ntnNij(ig1,nt12,nt22))
    bX(int1) = (g_ntnXij(ig1,nt11,nt21)+g_ntnXij(ig1,nt12,nt22))
    bSX(int1) = bS(int1) + bX(int1)
    bT(int1) = bN(int1) + bSX(int1)
    bSxu(int1) = bS(int1) * ubij(nt11,nt21)
    bNxu(int1) = bN(int1) * ubij(nt11,nt21)
    bXxu(int1) = bX(int1) * ubij(nt11,nt21)
    bSXxu(int1) = bSxu(int1) + bXxu(int1)
    bTxu(int1) = bNxu(int1) + bSXxu(int1)
    write(*,*) 'IAM\GEN: bS, bS x ubij   = ', bS(int1), bSxu(int1)
    write(*,*) 'IAM\GEN: bN, bN x ubij   = ', bN(int1), bNxu(int1)
    write(*,*) 'IAM\GEN: bX, bX x ubij   = ', bX(int1), bXxu(int1)
    write(*,*) 'IAM\GEN: bSX, bSX x ubij = ', bSX(int1), bSXxu(int1)
    write(*,*) 'IAM\GEN: bT, bT x ubij   = ', bT(int1), bTxu(int1)
    write(*,*) 'IAM\GEN: AT>TA'
    int1 = 3
    nt11 = 1
    nt21 = 2
    nt12 = 2
    nt22 = 1
    bS(int1) = (g_ntnSij(ig1,nt11,nt21)+g_ntnSij(ig1,nt12,nt22))
    bN(int1) = (g_ntnNij(ig1,nt11,nt21)+g_ntnNij(ig1,nt12,nt22))
    bX(int1) = (g_ntnXij(ig1,nt11,nt21)+g_ntnXij(ig1,nt12,nt22))
    bSX(int1) = bS(int1) + bX(int1)
    bT(int1) = bN(int1) + bSX(int1)
    bSxu(int1) = bS(int1) * ubij(nt11,nt21)
    bNxu(int1) = bN(int1) * ubij(nt11,nt21)
    bXxu(int1) = bX(int1) * ubij(nt11,nt21)
    bSXxu(int1) = bSxu(int1) + bXxu(int1)
    bTxu(int1) = bNxu(int1) + bSXxu(int1)
    write(*,*) 'IAM\GEN: bS, bS x ubij   = ', bS(int1), bSxu(int1)
    write(*,*) 'IAM\GEN: bN, bN x ubij   = ', bN(int1), bNxu(int1)
    write(*,*) 'IAM\GEN: bX, bX x ubij   = ', bX(int1), bXxu(int1)
    write(*,*) 'IAM\GEN: bSX, bSX x ubij = ', bSX(int1), bSXxu(int1)
    write(*,*) 'IAM\GEN: bT, bT x ubij   = ', bT(int1), bTxu(int1)
    write(*,*) 'IAM\GEN: GC>TA'
    int1 = 4
    nt11 = 4
    nt21 = 2
    nt12 = 3
    nt22 = 1
    bS(int1) = (g_ntnSij(ig1,nt11,nt21)+g_ntnSij(ig1,nt12,nt22))
    bN(int1) = (g_ntnNij(ig1,nt11,nt21)+g_ntnNij(ig1,nt12,nt22))
    bX(int1) = (g_ntnXij(ig1,nt11,nt21)+g_ntnXij(ig1,nt12,nt22))
    bSX(int1) = bS(int1) + bX(int1)
    bT(int1) = bN(int1) + bSX(int1)
    bSxu(int1) = bS(int1) * ubij(nt11,nt21)
    bNxu(int1) = bN(int1) * ubij(nt11,nt21)
    bXxu(int1) = bX(int1) * ubij(nt11,nt21)
    bSXxu(int1) = bSxu(int1) + bXxu(int1)
    bTxu(int1) = bNxu(int1) + bSXxu(int1)
    write(*,*) 'IAM\GEN: bS, bS x ubij   = ', bS(int1), bSxu(int1)
    write(*,*) 'IAM\GEN: bN, bN x ubij   = ', bN(int1), bNxu(int1)
    write(*,*) 'IAM\GEN: bX, bX x ubij   = ', bX(int1), bXxu(int1)
    write(*,*) 'IAM\GEN: bSX, bSX x ubij = ', bSX(int1), bSXxu(int1)
    write(*,*) 'IAM\GEN: bT, bT x ubij   = ', bT(int1), bTxu(int1)
    write(*,*) 'IAM\GEN: AT>CG'
    int1 = 5
    nt11 = 1
    nt21 = 3
    nt12 = 2
    nt22 = 4
    bS(int1) = (g_ntnSij(ig1,nt11,nt21)+g_ntnSij(ig1,nt12,nt22))
    bN(int1) = (g_ntnNij(ig1,nt11,nt21)+g_ntnNij(ig1,nt12,nt22))
    bX(int1) = (g_ntnXij(ig1,nt11,nt21)+g_ntnXij(ig1,nt12,nt22))
    bSX(int1) = bS(int1) + bX(int1)
    bT(int1) = bN(int1) + bSX(int1)
    bSxu(int1) = bS(int1) * ubij(nt11,nt21)
    bNxu(int1) = bN(int1) * ubij(nt11,nt21)
    bXxu(int1) = bX(int1) * ubij(nt11,nt21)
    bSXxu(int1) = bSxu(int1) + bXxu(int1)
    bTxu(int1) = bNxu(int1) + bSXxu(int1)
    write(*,*) 'IAM\GEN: bS, bS x ubij   = ', bS(int1), bSxu(int1)
    write(*,*) 'IAM\GEN: bN, bN x ubij   = ', bN(int1), bNxu(int1)
    write(*,*) 'IAM\GEN: bX, bX x ubij   = ', bX(int1), bXxu(int1)
    write(*,*) 'IAM\GEN: bSX, bSX x ubij = ', bSX(int1), bSXxu(int1)
    write(*,*) 'IAM\GEN: bT, bT x ubij   = ', bT(int1), bTxu(int1)
    write(*,*) 'IAM\GEN: GC>CG'
    int1 = 6
    nt11 = 4
    nt21 = 3
    nt12 = 3
    nt22 = 4
    bS(int1) = (g_ntnSij(ig1,nt11,nt21)+g_ntnSij(ig1,nt12,nt22))
    bN(int1) = (g_ntnNij(ig1,nt11,nt21)+g_ntnNij(ig1,nt12,nt22))
    bX(int1) = (g_ntnXij(ig1,nt11,nt21)+g_ntnXij(ig1,nt12,nt22))
    bSX(int1) = bS(int1) + bX(int1)
    bT(int1) = bN(int1) + bSX(int1)
    bSxu(int1) = bS(int1) * ubij(nt11,nt21)
    bNxu(int1) = bN(int1) * ubij(nt11,nt21)
    bXxu(int1) = bX(int1) * ubij(nt11,nt21)
    bSXxu(int1) = bSxu(int1) + bXxu(int1)
    bTxu(int1) = bNxu(int1) + bSXxu(int1)
    write(*,*) 'IAM\GEN: bS, bS x ubij   = ', bS(int1), bSxu(int1)
    write(*,*) 'IAM\GEN: bN, bN x ubij   = ', bN(int1), bNxu(int1)
    write(*,*) 'IAM\GEN: bX, bX x ubij   = ', bX(int1), bXxu(int1)
    write(*,*) 'IAM\GEN: bSX, bSX x ubij = ', bSX(int1), bSXxu(int1)
    write(*,*) 'IAM\GEN: bT, bT x ubij   = ', bT(int1), bTxu(int1)
    write(*,*) 'IAM\GEN: A+T>G+C'
    int1 = 1
    int2 = 5
    write(*,*) 'IAM\GEN: bS, bS x ubij   = ', bS(int1)+bS(int2), bSxu(int1)+bSxu(int2)
    write(*,*) 'IAM\GEN: bN, bN x ubij   = ', bN(int1)+bN(int2), bNxu(int1)+bNxu(int2)
    write(*,*) 'IAM\GEN: bX, bX x ubij   = ', bX(int1)+bX(int2), bXxu(int1)+bXxu(int2)
    write(*,*) 'IAM\GEN: bSX, bSX x ubij = ', bSX(int1)+bSX(int2), bSXxu(int1)+bSXxu(int2)
    write(*,*) 'IAM\GEN: bT, bT x ubij   = ', bT(int1)+bT(int2), bTxu(int1)+bTxu(int2)
    write(*,*) 'IAM\GEN: G+C>A+T'
    int1 = 2
    int2 = 4
    write(*,*) 'IAM\GEN: bS, bS x ubij   = ', bS(int1)+bS(int2), bSxu(int1)+bSxu(int2)
    write(*,*) 'IAM\GEN: bN, bN x ubij   = ', bN(int1)+bN(int2), bNxu(int1)+bNxu(int2)
    write(*,*) 'IAM\GEN: bX, bX x ubij   = ', bX(int1)+bX(int2), bXxu(int1)+bXxu(int2)
    write(*,*) 'IAM\GEN: bSX, bSX x ubij = ', bSX(int1)+bSX(int2), bSXxu(int1)+bSXxu(int2)
    write(*,*) 'IAM\GEN: bT, bT x ubij   = ', bT(int1)+bT(int2), bTxu(int1)+bTxu(int2)
    
    write(*,*) 'IAM\GEN: DNA CNP requirement:'
    write(*,*) 'IAM\GEN: q0DNAC =   ',  g_q0DNAC(ig1)
    write(*,*) 'IAM\GEN: q0DNAN =   ',  g_q0DNAN(ig1)
    write(*,*) 'IAM\GEN: q0DNAP =   ',  g_q0DNAP(ig1)
enddo ! ig1
write(*,*) 'IAM\GEN: uEg min, max', uEgmin, uEgmax

!
! --- process proteins ---
!
write(*,*) 'IAM\PRO:'
write(*,*) 'IAM\PRO: Calculating protein properties:'

do ig1 = 1, ngt

write(*,*) 'IAM\PRO:'
write(*,*) 'IAM\PRO: ig1 = ', ig1

write(*,*) 'IAM\PRO: Watch protein:'
write(*,*) 'IAM\PRO: start = ', g_pntstart(ig1,ipw)
write(*,*) 'IAM\PRO: stop = ', g_pntstop(ig1,ipw)
write(*,*) 'IAM\PRO: other = ', g_pother(ig1,ipw)

!
! adjust g_pgamma
!
pgammat = 0.
do ip1 = 1, g_np(ig1)
    pgammat = pgammat + g_pgamma(ig1,ip1)
enddo
pgammat = pgammat / dble(g_np(ig1))
!write(*,*) 'IAM\PRO: input ave(pgamma) = ', pgammat
do ip1 = 1, g_np(ig1)
    g_pgamma(ig1,ip1) = g_pgamma(ig1,ip1) / pgammat
enddo
!pgammat = 0.
!do ip1 = 1, g_np(ig1)
!    pgammat = pgammat + g_pgamma(ig1,ip1)
!enddo
!pgammat = pgammat / dble(g_np(ig1))
!write(*,*) 'IAM\PRO: adjusted ave(pgamma) = ', pgammat

!
! g_pother & min. quotas
!
potherf = 0.
g_q0aaC(ig1) = 0.
g_q0aaN(ig1) = 0.
q0aaCw = 0.
q0aaNw = 0.
g_naat(ig1) = 0.
naaxw = 0
naaz = 0
naazw = 0
ncodz = 0
g_naa(ig1,:) = 0
g_nntt(ig1,:) = 0.
g_nnttt(ig1) = 0
g_nnttw = 0.
g_nntttw = 0
pntnSt = 0.
pntnNt = 0.
do ip1 = 1, g_np(ig1)
    !write(*,*) 'IAM\PRO: ip1 = ', ip1
    !write(*,*) 'IAM\PRO: start, stop, other = ', g_pntstart(ig1,ip1), g_pntstop(ig1,ip1), g_pother(ig1,ip1)
    if (g_pother(ig1,ip1).eq.-1) then
        potherf = potherf + 1
    endif
    if (g_pntstop(ig1,ip1).gt.nnt) then
        write(*,*) "IAM\PRO: ERR: pntstop > nnt.", g_pntstop(ig1,ip1), nnt
        goto 9100
    endif
    if (g_pother(ig1,ip1).eq.-1) then
        pntstartx = g_pntstart(ig1,ip1)
        pntstopx = g_pntstop(ig1,ip1)
        pntstepx = 3
    else
        pntstartx = g_pntstop(ig1,ip1)
        pntstopx = g_pntstart(ig1,ip1)
        pntstepx = -3
    endif
    do int1 = pntstartx, pntstopx, pntstepx
        !write(*,*) 'IAM\PRO: int1 = ', int1
        if (g_pother(ig1,ip1).eq.-1) then
            nta(1) = g_nt(ig1,int1)
            nta(2) = g_nt(ig1,int1+1)
            nta(3) = g_nt(ig1,int1+2)
            cod1 = nta(1)*100+nta(2)*10+nta(3)*1
        else ! pother
            nta(1) = g_nt(ig1,int1)
            nta(2) = g_nt(ig1,int1-1)
            nta(3) = g_nt(ig1,int1-2)
            cod1 = int_opo(nta(1))*100+int_opo(nta(2))*10+int_opo(nta(3))*1
        endif ! pother
        iaa1 = gc(cod1)
        !write(*,*) 'IAM\PRO: cod1, iaa1 = ', cod1, iaa1
        if (iaa1.le.20) then
            g_nntt(ig1,g_nt(ig1,int1)) = g_nntt(ig1,g_nt(ig1,int1)) + g_pgamma(ig1,ip1)
            g_nntt(ig1,g_nt(ig1,int1+1)) = g_nntt(ig1,g_nt(ig1,int1+1)) + g_pgamma(ig1,ip1)
            g_nntt(ig1,g_nt(ig1,int1+2)) = g_nntt(ig1,g_nt(ig1,int1+2)) + g_pgamma(ig1,ip1)
            g_nnttt(ig1) = g_nnttt(ig1) + 3. * g_pgamma(ig1,ip1)
            g_naat(ig1) = g_naat(ig1) + 1
            naaz(iaa1) = naaz(iaa1) + 1
            g_naa(ig1,iaa1) = g_naa(ig1,iaa1) + 1
            g_q0aaC(ig1) = g_q0aaC(ig1) + g_pgamma(ig1,ip1) * raaC(iaa1)
            g_q0aaN(ig1) = g_q0aaN(ig1) + g_pgamma(ig1,ip1) * raaN(iaa1)
            if (ip1.eq.ipw) then
                g_nnttw(ig1,g_nt(ig1,int1)) = g_nnttw(ig1,g_nt(ig1,int1)) + g_pgamma(ig1,ip1)
                g_nnttw(ig1,g_nt(ig1,int1+1)) = g_nnttw(ig1,g_nt(ig1,int1+1)) + g_pgamma(ig1,ip1)
                g_nnttw(ig1,g_nt(ig1,int1+2)) = g_nnttw(ig1,g_nt(ig1,int1+2)) + g_pgamma(ig1,ip1)
                g_nntttw(ig1) = g_nntttw(ig1) + 3. * g_pgamma(ig1,ip1)
                naaxw = naaxw + 1
                naazw(iaa1) = naazw(iaa1) + 1
                q0aaCw = q0aaCw + g_pgamma(ig1,ip1) * raaC(iaa1)
                q0aaNw = q0aaNw + g_pgamma(ig1,ip1) * raaN(iaa1)
            endif
            ncodz(int_icod(cod1)) = ncodz(int_icod(cod1)) + 1
        endif
    enddo ! int1
    g_pntn(ig1,ip1) = g_pntnS(ig1,ip1) + g_pntnN(ig1,ip1)
    pntnt = pntnt + g_pntn(ig1,ip1)
    pntnSt = pntnSt + g_pntnS(ig1,ip1)
    pntnNt = pntnNt + g_pntnN(ig1,ip1)
enddo ! ip1

if (g_nnttt(ig1).eq.0) then
    write(*,*) 'IAM\PRO: ERR. g_nnttt = 0.'
    goto 9100
endif

write(*,*) 'IAM\PRO:'
write(*,*) 'IAM\PRO: Transcript pool composition.'
write(*,*) 'IAM\PRO: All proteins:'
write(*,*) 'IAM\PRO: n, f A =   ', g_nntt(ig1,1), g_nntt(ig1,1) / g_nnttt(ig1)
write(*,*) 'IAM\PRO: n, f T =   ', g_nntt(ig1,2), g_nntt(ig1,2) / g_nnttt(ig1)
write(*,*) 'IAM\PRO: n, f C =   ', g_nntt(ig1,3), g_nntt(ig1,3) / g_nnttt(ig1)
write(*,*) 'IAM\PRO: n, f G =   ', g_nntt(ig1,4), g_nntt(ig1,4) / g_nnttt(ig1)
write(*,*) 'IAM\PRO: n, f GC =  ', (g_nntt(ig1,3)+g_nntt(ig1,4)), (g_nntt(ig1,3)+g_nntt(ig1,4)) / g_nnttt(ig1)
write(*,*) 'IAM\PRO: Watch protein: ipw = ', ipw
write(*,*) 'IAM\PROw: n, f A =   ', g_nnttw(ig1,1), g_nnttw(ig1,1) / g_nntttw(ig1)
write(*,*) 'IAM\PROw: n, f T =   ', g_nnttw(ig1,2), g_nnttw(ig1,2) / g_nntttw(ig1)
write(*,*) 'IAM\PROw: n, f C =   ', g_nnttw(ig1,3), g_nnttw(ig1,3) / g_nntttw(ig1)
write(*,*) 'IAM\PROw: n, f G =   ', g_nnttw(ig1,4), g_nnttw(ig1,4) / g_nntttw(ig1)
write(*,*) 'IAM\PROw: n, f GC =  ', (g_nnttw(ig1,3)+g_nnttw(ig1,4)), (g_nnttw(ig1,3)+g_nnttw(ig1,4)) / g_nntttw(ig1)

write(*,*) 'IAM\PRO:'
! note: fGCt and fGCt0 are temporary variables, used only here
fGCt = (g_nntt(ig1,3)+g_nntt(ig1,4)) / g_nnttt(ig1)
if (oLothIC.eq.1) then
    write(*,*) 'IAM\PRO: fGCt0 based on genome.'
    fGCt0 = fGCt
else
    write(*,*) 'IAM\PRO: fGCt0 based on input.'
endif
Loth = 1.
if (oLoth.eq.1) then
    Loth = 1. + (fGCt - fGCt0) * sGCt
    write(*,*) 'IAM\PRO: Loth = ', Loth
endif
g_Vmax0C(ig1) = Vmax0C * Loth
g_Vmax0N(ig1) = Vmax0N * Loth
g_Vmax0P(ig1) = Vmax0P * Loth

potherf = potherf / dble(g_np(ig1))
write(*,*) 'IAM\PRO:'
write(*,*) 'IAM\PRO: other fraction = ', potherf
g_q0aaC(ig1) = m0aa / dble(g_naat(ig1)) * g_q0aaC(ig1)
g_q0aaN(ig1) = m0aa / dble(g_naat(ig1)) * g_q0aaN(ig1)
q0aaCw = m0aa / dble(naaxw) * q0aaCw
q0aaNw = m0aa / dble(naaxw) * q0aaNw
write(*,*) 'IAM\PRO: np = ', g_np(ig1)
write(*,*) 'IAM\PRO: Sites:'
write(*,*) 'IAM\PRO: Average across proteins:'
write(*,*) 'IAM\PRO: T = ', pntnt/g_np(ig1)
write(*,*) 'IAM\PRO: S = ', pntnSt/g_np(ig1)
write(*,*) 'IAM\PRO: N = ', pntnNt/g_np(ig1)
write(*,*) 'IAM\PRO: Watch protein:'
write(*,*) 'IAM\PRO: T = ', g_pntn(ig1,ipw)
write(*,*) 'IAM\PRO: S = ', g_pntnS(ig1,ipw)
write(*,*) 'IAM\PRO: N = ', g_pntnN(ig1,ipw)
write(*,*) 'IAM\PRO: naat = ', g_naat(ig1)
write(*,*) 'IAM\PRO: naat*3/nnt = ', g_naat(ig1)*3./dble(nnt)
write(*,*) 'IAM\PRO: q0aaC = ', g_q0aaC(ig1)
write(*,*) 'IAM\PRO: q0aaN = ', g_q0aaN(ig1)
write(*,*) 'IAM\PRO: q0aaCw = ', q0aaCw
write(*,*) 'IAM\PRO: q0aaNw = ', q0aaNw
write(*,*) 'IAM\PRO: Amino acid composition:'
do iaa1 = 1, 20
    write(*,*) 'IAM\PRO: ', iaa1, int_saa(iaa1), naaz(iaa1), dble(naaz(iaa1))/dble(g_naat(ig1))*100.
enddo
do iaa1 = 1, 20
    write(*,*) 'IAM\PROw: ', iaa1, int_saa(iaa1), naazw(iaa1), dble(naazw(iaa1))/dble(naaxw)*100.
enddo
write(*,*) 'IAM\PRO: Codon composition:'
do icod1 = 1, 64
    write(*,*) 'IAM\PRO: ', icod1, int_scod(icod1), ncodz(icod1), dble(ncodz(icod1))/dble(g_naat(ig1))*100.
enddo

write(*,*) 'IAM\INIT: q0C = ', q0othC + g_q0DNAC(ig1) + g_q0aaC(ig1)
write(*,*) 'IAM\INIT: q0N = ', q0othN + g_q0DNAN(ig1) + g_q0aaN(ig1)
write(*,*) 'IAM\INIT: q0P = ', q0othP + g_q0DNAP(ig1)

!
! spacing
!
pminspacing = dimnnt
pmaxspacing = - dimnnt
pminlength = dimnnt
pmaxlength = - dimnnt
do ip1 = 2, g_np(ig1)
    int2 = g_pntstart(ig1,ip1) - g_pntstop(ig1,ip1-1)
    if (int2.lt.pminspacing) pminspacing = int2
    if (int2.gt.pmaxspacing) pmaxspacing = int2
    int2 = g_pntstop(ig1,ip1) - g_pntstart(ig1,ip1)
    if (int2.lt.pminlength) pminlength = int2
    if (int2.gt.pmaxlength) pmaxlength = int2
enddo ! ip1
write(*,*) 'IAM\PRO: spacing min, max = ', pminspacing, pmaxspacing
write(*,*) 'IAM\PRO: length min, max = ', pminlength, pmaxlength

enddo ! ig1


!
! -------------------
! -------------------
! --- initial pop ---
! -------------------
! -------------------
!
write(*,*) 'IAM\INIT:'
write(*,*) 'IAM\INIT: Initial population.'

nat = 0
SRT = 0.

if ((oHot.eq.0).or.(oHot.eq.1)) then

write(*,*) 'IAM\INIT: Cold start.'
    
nic_sp = dble(nic) / dble(dimnsp)

do isp1 = 1, dimnsp
nat_sp(1,isp1) = 0
naf_sp(isp1) = 0
Popic_sp(isp1) = Popic / dble(dimnsp)
gsn_sp(isp1) = gsn
do ia1 = 1, nic_sp
    if ((nat_sp(1,isp1)+1).gt.dimna_sp) then
        write(*,*) 'IAM\INIT: ERR: nat_sp > dimna_sp.', nat_sp(1,isp1), dimna_sp
        write(*,*) 'IAM\INIT: nic_sp = ', nic_sp
        goto 9100
    endif
    nat_sp(1,isp1) = nat_sp(1,isp1) + 1
    nat = nat + 1

    a_sa(1,isp1,ia1) = 1
    a_sr(1,isp1,ia1) = Popic / dble(nic)
    SRT = SRT + a_sr(1,isp1,ia1)
    a_srm(1,isp1,ia1) = 0.
    a_srr(1,isp1,ia1) = 0.
    a_iar(1,isp1,ia1) = -9
    if (oStartDNA.eq.1) then
        ig1 = ia1
    else
        ig1 = 1
    endif
    a_gid(1,isp1,ia1) = ig1
    a_cn(1,isp1,ia1) = 0
    !do ic1 = 1, dimnc
    !    a_ck(ic1,isp1,ia1,1) = 0
    !    a_cnt(ic1,isp1,ia1,1) = 0
    !enddo
    do nt1 = 1, 4
        a_ntnk(1,isp1,ia1,nt1) = g_ntnk(ig1,nt1)
    enddo
    do iaa1 = 1, 20
        a_naa(1,isp1,ia1,iaa1) = g_naa(ig1,iaa1)
        do iaa2 = 1, 20
            a_naax(1,isp1,ia1,iaa1,iaa2) = 0
        enddo
    enddo
    a_uEg(1,isp1,ia1) = g_uEg(ig1)
    a_q0DNAC(1,isp1,ia1) = g_q0DNAC(ig1)
    a_q0DNAN(1,isp1,ia1) = g_q0DNAN(ig1)
    a_q0DNAP(1,isp1,ia1) = g_q0DNAP(ig1)
    a_q0aaC(1,isp1,ia1) = g_q0aaC(ig1)
    a_q0aaN(1,isp1,ia1) = g_q0aaN(ig1)
    a_VmaxC(1,isp1,ia1) = g_Vmax0C(ig1)
    a_VmaxN(1,isp1,ia1) = g_Vmax0N(ig1)
    a_VmaxP(1,isp1,ia1) = g_Vmax0P(ig1)
    if (oSlim.eq.0) then
        call new_sn(1,isp1,ia1)
        a_rsn(1,isp1,ia1) = a_sn(1,isp1,ia1)
        a_nd(1,isp1,ia1) = 0
        do nt1 = 1, 4
            a_ntnSk(1,isp1,ia1,nt1) = g_ntnSk(ig1,nt1)
            a_ntnNk(1,isp1,ia1,nt1) = g_ntnNk(ig1,nt1)
            a_ntnXk(1,isp1,ia1,nt1) = g_ntnXk(ig1,nt1)
        enddo
        a_nmt(1,isp1,ia1) = 0
        a_nmS(1,isp1,ia1) = 0
        a_nmN(1,isp1,ia1) = 0
        a_nrt(1,isp1,ia1) = 0
        a_nrS(1,isp1,ia1) = 0
        a_nrN(1,isp1,ia1) = 0
        do ic1 = 1, dimnc
            a_cXSN(ic1,isp1,ia1,1) = 0
            a_cmr(ic1,isp1,ia1,1) = 0
        enddo
        do nt1 = 1, 4
            do nt2 = 1, 4
                a_nm(1,isp1,ia1,nt1,nt2) = 0
                a_nr(1,isp1,ia1,nt1,nt2) = 0
            enddo
        enddo
        do nt1 = 1, 4
            a_nntt(1,isp1,ia1,nt1) = g_nntt(ig1,nt1)
        enddo
        a_VC(1,isp1,ia1) = g_Vmax0C(ig1)
        a_VN(1,isp1,ia1) = g_Vmax0N(ig1)
        a_VP(1,isp1,ia1) = g_Vmax0P(ig1)
        q0C = q0othC + g_q0DNAC(ig1) + g_q0aaC(ig1)
        q0N = q0othN + g_q0DNAN(ig1) + g_q0aaN(ig1)
        q0P = q0othP + g_q0DNAP(ig1)
        if (oHANS.eq.0) then
            r = r4_uni(idum_sp(isp1))
            a_qC(1,isp1,ia1) = q0C + r * q0C
            a_qN(1,isp1,ia1) = q0N + r * q0N
            a_qP(1,isp1,ia1) = q0P + r * q0P
        else
            qaveC = q0C * log(2.) * 2.
            qaveN = q0N * log(2.) * 2.
            qaveP = q0P * log(2.) * 2.
            a_qC(1,isp1,ia1) = qaveC
            a_qN(1,isp1,ia1) = qaveN
            a_qP(1,isp1,ia1) = qaveP
        endif
        a_CLimit(1,isp1,ia1) = 2.
        a_NLimit(1,isp1,ia1) = 2.
        a_PLimit(1,isp1,ia1) = 2.
        a_tb(1,isp1,ia1) = 0.
        a_tg(1,isp1,ia1) = log(2.) / (g_Vmax0C(ig1) / a_qC(1,isp1,ia1))
    endif ! oSlim
        
enddo ! ia1

if (gsn_sp(isp1).gt.gsn) then
    gsn = gsn_sp(isp1)
endif

enddo ! isp1

else ! oHot
    
write(*,*) 'IAM\INIT: Hot start.'

do isp1 = 1, dimnsp
    nat_sp(1,isp1) = 0
    naf_sp(isp1) = 0
    SRT_sp(isp1) = 0.
    Popic_sp(isp1) = Popic / dble(dimnsp)
    gsn_sp(isp1) = gsn
enddo

natHOT = 0
SRTHOT = 0.
ncomp = 100
do icomp = 1, ncomp
    nat_comp(icomp) = 0
    SRT_comp(icomp) = 0.
    do isp1 = 1, dimnsp
        nat_comp_sp(icomp,isp1) = 0
        SRT_comp_sp(icomp,isp1) = 0.
    enddo
enddo
ncomp = 0

open(421,file='IAM_H_POP.bin', form='unformatted')
do while(not(eof(421)))
    read(421) &
        isp1, &
        ia1, & !< different from write
        a_sa(1,isp1,ia1), &
        a_sr(1,isp1,ia1), &
        a_srm(1,isp1,ia1), & 
        a_srr(1,isp1,ia1), &
        a_iar(1,isp1,ia1), &
        a_gid(1,isp1,ia1), &
        a_cn(1,isp1,ia1)
    if (ia1.gt.dimna_sp) then
        write(*,*) 'IAM\HOT: ERR: ia1 > dimna_sp.', ia1, dimna_sp
        goto 9300
    endif
    if (a_sa(1,isp1,ia1).lt.1) then
        write(*,*) 'IAM\HOT: ERR. sa<1 read.'
        write(*,*) 'IAM\HOT: isp1, ia1 = ', isp1, ia1
        write(*,*) 'IAM\HOT: sa = ', a_sa(1,isp1,ia1)
    endif
    natHOT = natHOT + 1
    SRTHOT = SRTHOT + a_sr(1,isp1,ia1)
    read(421) (a_ck(ic1,isp1,ia1,1), ic1 = 1, a_cn(1,isp1,ia1))
    read(421) (a_cnt(ic1,isp1,ia1,1), ic1 = 1, a_cn(1,isp1,ia1))
    read(421) (a_ntnk(1,isp1,ia1,nt1), nt1 = 1, 4)
    do ic1 = 1, a_cn(1,isp1,ia1)
        if ((a_cnt(ic1,isp1,ia1,1).lt.1).or.(a_cnt(ic1,isp1,ia1,1).gt.4)) then
            write(*,*) 'IAM\HOT: ERR. Invalid nt read.'
            write(*,*) 'IAM\HOT: isp1, ia1, ic1 = ', isp1, ia1, ic1
            write(*,*) 'IAM\HOT: nt = ', a_cnt(ic1,isp1,ia1,1)
            goto 9100
        endif
    enddo
    read(421) (a_naa(1,isp1,ia1,iaa1), iaa1 = 1, 20)
    do iaa1 = 1, 20
        read(421) (a_naax(1,isp1,ia1,iaa1,iaa2), iaa2 = 1, 20)
    enddo
    read(421) &
        a_uEg(1, isp1,ia1), &
        a_q0DNAC(1,isp1,ia1), &
        a_q0DNAN(1,isp1,ia1), &
        a_q0DNAP(1,isp1,ia1), &
        a_q0aaC(1,isp1,ia1), &
        a_q0aaN(1,isp1,ia1), &
        a_VmaxC(1,isp1,ia1), &
        a_VmaxN(1,isp1,ia1), &
        a_VmaxP(1,isp1,ia1)
    if ( &
        isnan(a_uEg(1, isp1,ia1)).or. &
        isnan(a_q0DNAC(1,isp1,ia1)).or. &
        isnan(a_q0DNAN(1,isp1,ia1)).or. &
        isnan(a_q0DNAP(1,isp1,ia1)).or. &
        isnan(a_q0aaC(1,isp1,ia1)).or. &
        isnan(a_q0aaN(1,isp1,ia1)).or. &
        isnan(a_VmaxC(1,isp1,ia1)).or. &
        isnan(a_VmaxN(1,isp1,ia1)).or. &
        isnan(a_VmaxP(1,isp1,ia1))) then
        write(*,*) 'IAM\HOT: ERR. nan read.'
        write(*,*) 'IAM\HOT: isp1, ia1 = ', isp1, ia1
        write(*,*) 'IAM\HOT: ', &
            a_uEg(1, isp1,ia1), &
            a_q0DNAC(1,isp1,ia1), &
            a_q0DNAN(1,isp1,ia1), &
            a_q0DNAP(1,isp1,ia1), &
            a_q0aaC(1,isp1,ia1), &
            a_q0aaN(1,isp1,ia1), &
            a_VmaxC(1,isp1,ia1), &
            a_VmaxN(1,isp1,ia1), &
            a_VmaxP(1,isp1,ia1)
        endif
if (oSlim.eq.0) then
        read(421) &
            a_sn(1,isp1,ia1), &
            a_rsn(1,isp1,ia1), &
            a_nd(1,isp1,ia1), &
            (a_ntnSk(1,isp1,ia1,nt1), nt1 = 1, 4), &
            (a_ntnNk(1,isp1,ia1,nt1), nt1 = 1, 4), &
            (a_ntnXk(1,isp1,ia1,nt1), nt1 = 1, 4), &
            a_nmt(1,isp1,ia1), &
            a_nmS(1,isp1,ia1), &
            a_nmN(1,isp1,ia1), &
            a_nrt(1,isp1,ia1), &
            a_nrS(1,isp1,ia1), &
            a_nrN(1,isp1,ia1)
        ! check
        icomp = a_srr(1,isp1,ia1)
        if (icomp.gt.ncomp) then
            ncomp = icomp
        endif
        if (icomp.ge.1) then
            nat_comp(icomp) = nat_comp(icomp) + 1
            SRT_comp(icomp) = SRT_comp(icomp) + a_sr(1,isp1,ia1)
            nat_comp_sp(icomp,isp1) = nat_comp_sp(icomp,isp1) + 1
            SRT_comp_sp(icomp,isp1) = SRT_comp_sp(icomp,isp1) + a_sr(1,isp1,ia1)
        endif
        if (a_srm(1,isp1,ia1).lt.0.) then
            write(*,*) 'IAM\HOT\CHK: in icomp from, SN = ', nint(a_srr(1,isp1,ia1)), a_sn(1,isp1,ia1)
        endif
        a_srm(1,isp1,ia1) = 0.
        a_srr(1,isp1,ia1) = 0.
        read(421) (a_cXSN(ic1,isp1,ia1,1), ic1 = 1, a_cn(1,isp1,ia1))
        read(421) (a_cmr(ic1,isp1,ia1,1), ic1 = 1, a_cn(1,isp1,ia1))
        do ic1 = 1, a_cn(1,isp1,ia1)
            if ((a_cXSN(ic1,isp1,ia1,1).lt.0).or.(a_cmr(ic1,isp1,ia1,1).lt.0)) then
                write(*,*) 'IAM\HOT: ERR. Invalid cXSN or cmr read.'
                write(*,*) 'IAM\HOT: isp1, ia1, ic1 = ', isp1, ia1, ic1
                write(*,*) 'IAM\HOT: cXSN = ', a_cXSN(ic1,isp1,ia1,1)
                write(*,*) 'IAM\HOT: cmr = ', a_cmr(ic1,isp1,ia1,1)
                goto 9100
            endif
        enddo
        read(421) ((a_nm(1,isp1,ia1,nt1,nt2), nt1 = 1, 4), nt2 = 1, 4)
        read(421) ((a_nr(1,isp1,ia1,nt1,nt2), nt1 = 1, 4), nt2 = 1, 4)
        read(421) (a_nntt(1,isp1,ia1,nt1), nt1 = 1, 4)
        read(421)  &
            a_VC(1,isp1,ia1), &
            a_VN(1,isp1,ia1), &
            a_VP(1,isp1,ia1), &
            a_qC(1,isp1,ia1), &
            a_qN(1,isp1,ia1), &
            a_qP(1,isp1,ia1), &
            a_CLimit(1,isp1,ia1), &
            a_NLimit(1,isp1,ia1), &
            a_PLimit(1,isp1,ia1), &
            a_tb(1,isp1,ia1), &
            a_tg(1,isp1,ia1)
    endif

    nat_sp(1,isp1) = nat_sp(1,isp1) + 1
    SRT_sp(isp1) = SRT_sp(isp1) + a_sr(1,isp1,ia1)
    nat = nat + 1
    if (a_sn(1,isp1,ia1).gt.gsn_sp(isp1)) then
        gsn_sp(isp1) = a_sn(1,isp1,ia1)
    endif
    
enddo ! eof
close(421)

write(*,*) 'IAM\INIT: ncomp = ', ncomp
write(*,*) 'IAM\INIT: natHOT, SRTHOT = ', natHOT, SRTHOT
write(*,*) 'IAM\INIT: icomp, nat, %:'
do icomp = 1, ncomp
    write(*,*) 'IAM\INIT:', icomp, nat_comp(icomp), dble(nat_comp(icomp))/dble(natHOT)*100.
enddo
write(*,*) 'IAM\INIT: icomp, SRT, %:'
do icomp = 1, ncomp
    write(*,*) 'IAM\INIT:', icomp, SRT_comp(icomp), SRT_comp(icomp)/SRTHOT*100.
enddo

write(*,*) 'IAM\INIT: isp1, nat_sp, SRT_sp:'
do isp1 = 1, dimnsp
    write(*,*) 'IAM\INIT: ', isp1, nat_sp(1,isp1), SRT_sp(isp1)
    if (gsn_sp(isp1).gt.gsn) then
        gsn = gsn_sp(isp1)
    endif
enddo

write(*,*) 'IAM\INIT: SRTHOT, Popic = ', SRTHOT, Popic
if (oHots.eq.1) then
write(*,*) 'IAM\INIT: Adjusting SR for Popic.'
ia2 = 0
do isp1 = 1, dimnsp
    do ia1 = 1, nat_sp(1,isp1)
        if (a_sa(1,isp1,ia1).ge.1) then
            a_sr(1,isp1,ia1) = a_sr(1,isp1,ia1) * Popic / SRTHOT
            if (a_sr(1,isp1,ia1).lt.1.) then
                 call kill_agent(isp1,ia1)
                 ia2 = ia2 + 1
            endif
        endif
    enddo
enddo
write(*,*) 'IAM\INIT: Agents killed (n, %) = ', ia2, dble(ia2)/natHOT*100.
endif

!write(*,*) 'IAM\INIT: icomp, nat_comp_sp, SRT_comp_sp:'
!do isp1 = 1, dimnsp
!    write(*,*) 'IAM\INIT: isp1 = ', isp1
!    do icomp = 1, ncomp
!        write(*,*) 'IAM\INIT: ', icomp, nat_comp_sp(icomp,isp1), SRT_comp_sp(icomp,isp1)
!    enddo
!enddo

endif ! oHot

naa = nat

write(*,*) 'IAM\INIT: Initial population done.'
write(*,*) 'IAM\INIT: nat = ', nat

do isp1 = 1, dimnsp
    VOL_sp(isp1) = VOL / dble(dimnsp)
enddo

!
! ---------------------
! --- set up output ---
! ---------------------
!
! --- cluster output file name string ---
!
if (cTime.gt.0) then
    write(*,*) 'IAM\CLU: '
    write(*,*) 'IAM\CLU: Setting up cluster output.'
    cString = '_t'//sap(cTime)//'_c'//sap(cComp)//'_i'//sap(cIter)
    write(*,*) 'IAM\CLU: cString = ', cString
endif

!
! --- open STA output file ---
!
if (cTime.gt.0) then
    open(241,file='IAM_O_STA'//cString//'.txt')
    open(242,file='IAM_O_STB'//cString//'.txt')
else
    open(241,file='IAM_O_STA.txt')
    open(242,file='IAM_O_STB.txt')
endif

tNoteNext = dtNote

tpnSTA = tpsSTA
ipSTAf = 0
ipSTAn = 0

nadt = 0.
nsrdt = 0.
fmutt = 0.
frect = 0.
fmutrecn = 0.
namt = 0.
nsrmt = 0.
nmTt = 0.
nmSt = 0.
nmNt = 0.
nmXt = 0.
do nt1 = 1, 4
    do nt2 = 1, 4
        nmTij(nt1,nt2) = 0.
        nmSij(nt1,nt2) = 0.
        nmNij(nt1,nt2) = 0.
        nmXij(nt1,nt2) = 0.
    enddo
enddo
nart = 0.
nsrrt = 0.
nrTt = 0.
nrSt = 0.
nrNt = 0.
nrXt = 0.
do nt1 = 1, 4
    do nt2 = 1, 4
        nrTij(nt1,nt2) = 0.
        nrSij(nt1,nt2) = 0.
        nrNij(nt1,nt2) = 0.
        nrXij(nt1,nt2) = 0.
    enddo
enddo
rdeltat = 0.
rnut = 0.
rn = 0.
do int1 = 1, 9
    rtype(int1) = 0.
enddo
rtypet = 0.
fpt = 0.
fpn2 = 0.
nawt = 0.
nsrwt = 0.

t_first = -9.
t_last = -9.
GC_first = -9.
GCT_first = -9.
GCS_first = -9.
GCN_first = -9.
GCX_first = -9.
GCXS_first = -9.
GC_last = -9.
GCT_last = -9.
GCS_last = -9.
GCN_last = -9.
GCX_last = -9.
GCXS_last = -9.
piS_first = -9.
piS_last = -9.
piN_first = -9.
piN_last = -9.
piM_first = -9.
piM_last = -9.
piR_first = -9.
piR_last = -9.
piT_first = -9.
piT_last = -9.
dNdS_first = -9.
dNdS_last = -9.
rm_first = -9.
rm_last = -9.
q0DNAC_first = -9.
q0DNAC_last = -9.
q0DNAN_first = -9.
q0DNAN_last = -9.
q0DNAP_first = -9.
q0DNAP_last = -9.
q0aaC_first = -9.
q0aaC_last = -9.
q0aaN_first = -9.
q0aaN_last = -9.
q0C_first = -9.
q0C_last = -9.
q0N_first = -9.
q0N_last = -9.
q0P_first = -9.
q0P_last = -9.
kg_first = -9.
kg_last = -9.

tpnGLO = tpsGLO
ipGLO = 0

ipCLK = 0

tpnPOP = tpsPOP
ipPOP = 0

tpnORF = tpsORF

!
! -----------------
! -----------------
! --- iteration ---
! -----------------
! -----------------
!
if (oIter.eq.1) then

    write(*,*) 'IAM\ITER: Saving population for iter.'
    do isp1 = 1, dimnsp
        nat_sp(2,isp1) = nat_sp(1,isp1) 
        do ia1 = 1, nat_sp(1,isp1)
            call copy_agent(1,isp1,ia1,2,isp1,ia1)
        enddo
    enddo

    iIter = 1

endif

cp_init_delta = cp_init_delta + omp_get_wtime() - cp_init_t1

8100 continue ! iteration start

if (oIter.eq.1) then
if (iIter.gt.1) then

    cp_iter_t1 = omp_get_wtime()
    
    write(*,*) 'IAM\ITER: Reverting to saved population.'
    !$OMP parallel num_threads(ntd) default(shared) PRIVATE(isp1, ia1)
    !$OMP DO
    do isp1 = 1, dimnsp
        nat_sp(1,isp1) = nat_sp(2,isp1) 
        naf_sp(isp1) = 0
        do ia1 = 1, nat_sp(2,isp1)
            call copy_agent(2,isp1,ia1,1,isp1,ia1)
        enddo
    enddo
    !$OMP END DO
    !$OMP END PARALLEL
    
    tpnSTA = tpsSTA + dtpSTA
    ipSTAf = 0

    cp_iter_delta = cp_iter_delta + omp_get_wtime() - cp_iter_t1

endif
endif

!
! ------------------
! ------------------
! --- time setup ---
! ------------------
! ------------------
!
t = 0.
th = 0.
iHOT = 0
tHOT = 0.
natHOT = 0
SRTHOT = 0.
if ((oHot.eq.2).or.(oHot.eq.3)) then
    open(411,file='IAM_H_OTH.txt')
    read(411,*) iHOT, tHOT, natHOT, SRTHOT
    close(411)
    write(*,*) 'IAM\HOT: iHOT, tHOT = ', iHOT, tHOT
    write(*,*) 'IAM\HOT: natHOT, SRTHOT = ', natHOT, SRTHOT
    iHOT = iHOT + 1
    th = tHOT
endif

ntx = max(1,nint(dtmix/dt))

!
! --- error handling ---
!
do isp1 = 1, dimnsp
    err_sp(isp1) = 0
enddo


!
! -----------------
! -----------------
! --- main loop ---
! -----------------
! -----------------
!
!
! -----------------------
! --- outer time loop ---
! -----------------------
!
do while (.true.)
!write(*,*) 'IAM: Outer time start: t, th = ', t, th

!
! -----------------------------
! --- update temporal input ---
! -----------------------------
!
if (t.ge.tt(it)) then

write(*,*) 'IAM\TEMP: Temporal input update:'
write(*,*) 'IAM\TEMP: it, t = ', it, t

5110 continue

fQ = tfQ(it)
fSinC = tfSinC(it)
fSinN = tfSinN(it)
fSinP = tfSinP(it)
fSCfix = tfSCfix(it)
fSNfix = tfSNfix(it)
fSPfix = tfSPfix(it)
fX1 = tfX1(it)
fX2 = tfX2(it)

!write(*,*) 'IAM\TEMP: fQ = ', fQ
!write(*,*) 'IAM\TEMP: fSinC = ', fSinC
!write(*,*) 'IAM\TEMP: fSinN = ', fSinN
!write(*,*) 'IAM\TEMP: fSinP = ', fSinP
!write(*,*) 'IAM\TEMP: fSCfix = ', fSCfix
!write(*,*) 'IAM\TEMP: fSNfix = ', fSNfix
!write(*,*) 'IAM\TEMP: fSPfix = ', fSPfix
!write(*,*) 'IAM\TEMP: fX1 = ', fX1
!write(*,*) 'IAM\TEMP: fX2 = ', fX2

Q = QB * fQ
do isp1 = 1, dimnsp
    Q_sp(isp1) = Q / dble(dimnsp)
enddo

SinC = SinCB * fSinC
SinN = SinNB * fSinN
SinP = SinPB * fSinP

if (oFixExt.eq.2) then
    do isp1 = 1, dimnsp
        SC_sp(isp1) = SC0 * fSCfix
        SN_sp(isp1) = SN0 * fSNfix
        SP_sp(isp1) = SP0 * fSPfix
    enddo
endif

if (oHANS.eq.1) then
if (oHANSs.eq.3) then
    isarkar = 0
endif
endif

it = it + 1

if (t.ge.tt(it)) then
    !write(*,*) 'IAM\TEMP: Skipping ahead...'
    goto 5110
endif
    
endif

!
! -----------------------------------------
! --- calc population size & limitation ---
! -----------------------------------------
!
if (oLP.eq.1) then
    P=0
    do isp1 = 1, dimnsp
    do ia1 = 1, nat_sp(1,isp1)
        if (a_sa(1,isp1,ia1).ge.1) then
            P = P + 1
        endif
    enddo
    enddo
    LP = max(1. - P / K, 0.)
else
    LP = 1.
endif

!
! ----------------------------------
! --- update fixed changes count ---
! ----------------------------------
!
5210 continue
jcg = 0
if (ncg_sp_mx.ge.dimnc) then
    write(*,*) 'IAM\FIX: Updating fixed changes count (ncg_sp_mx > dimnc).'
    jcg = 1
    icg = 1
    goto 5215
endif
if ((oFixOut.eq.1).and.(ipSTAf.eq.1)) then
    write(*,*) 'IAM\FIX: Updating fixed changes count (final).'
    jcg = 1
    goto 5215
endif
if ((oFixOut.eq.1).and.(t.ge.tpnSTA).and.(ipD.eq.1)) then
    write(*,*) 'IAM\FIX: Updating fixed changes count (output).'
    jcg = 1
    goto 5215
endif
5215 continue
if (oIter.eq.1) then
    write(*,*) 'IAM\ITER: Skipping fixed changes count.'
    goto 5220
endif
if (jcg.eq.1) then
    call fix_changes()
endif
5220 continue

!
! ----------------------
! --- update genomes ---
! ----------------------
!
if (icg.eq.1) then
cp_genomes_t1 = omp_get_wtime()
icg = 0
if (g_ncft.gt.0) then
    write(*,*) 'IAM\GEN: Updating genomes. t = ', t
    write(*,*) 'IAM\GEN: g_ncft = ', g_ncft
    do ig1 = 1, ngt
        ifx(ig1) = 0
    enddo
    do isp1 = 1, dimnsp
    do ia1 = 1, nat_sp(1,isp1)
        if (a_sa(1,isp1,ia1).ge.1) then
            ig1 = a_gid(1,isp1,ia1)
            if (g_ncf(ig1).eq.0) then
                goto 7202
            endif
            ic3 = 0
            if (a_cn(1,isp1,ia1).eq.fx_cn(ig1)) then
                a_cn(1,isp1,ia1) = 0
            else
                do ic1 = fx_cn(ig1) + 1, a_cn(1,isp1,ia1)
                    ic3 = ic3 + 1
                    a_ck(ic3,isp1,ia1,1) = a_ck(ic1,isp1,ia1,1)
                    a_cnt(ic3,isp1,ia1,1) = a_cnt(ic1,isp1,ia1,1)
                    a_cXSN(ic3,isp1,ia1,1) = a_cXSN(ic1,isp1,ia1,1)
                    a_cmr(ic3,isp1,ia1,1) = a_cmr(ic1,isp1,ia1,1)
                enddo
                a_cn(1,isp1,ia1) = a_cn(1,isp1,ia1) - fx_cn(ig1)
            endif
7202 continue
        endif ! sa
    enddo ! ia1
    enddo ! isp1
    do ig1 = 1, ngt
    do ic1 = 1, fx_cn(ig1)
         g_nt(ig1,fx_ck(ic1,ig1)) = fx_cnt(ic1,ig1)
         g_cmr(ig1,fx_ck(ic1,ig1)) = fx_cmr(ic1,ig1)
    enddo
    g_ncff(ig1) = g_ncff(ig1) + g_ncf(ig1)
    g_ncSff(ig1) = g_ncSff(ig1) + g_ncSf(ig1)
    g_ncNff(ig1) = g_ncNff(ig1) + g_ncNf(ig1)
    g_ncmff(ig1) = g_ncmff(ig1) + g_ncmf(ig1)
    g_ncrff(ig1) = g_ncrff(ig1) + g_ncrf(ig1)
    g_ncf(ig1) = 0
    g_ncSf(ig1) = 0
    g_ncNf(ig1) = 0
    g_ncmf(ig1) = 0
    g_ncrf(ig1) = 0
    enddo
endif    
cp_genomes_delta = cp_genomes_delta + omp_get_wtime() - cp_genomes_t1
endif

cp_output_t1 = omp_get_wtime()

!
! ------------------
! --- output STA ---
! ------------------
!
if (((t.ge.tpnSTA).and.(ipD.eq.1)).or.(ipSTAf.eq.1)) then
write(*,*) 'IAM\PSTA: t, tpnSTA = ', t, tpnSTA
ipSTAn = ipSTAn + 1

!
! --- statistics ---
!
! species, based on rsn
!
if (opSTA.eq.1) then
    ns = 0
    do isp1 = 1, dimnsp
    do ia1 = 1, nat_sp(1,isp1)
        if (a_sa(1,isp1,ia1).ge.1) then
            do is1 = 1, ns
                if (a_rsn(1,isp1,ia1).eq.srsn(is1)) then !this guy is not unique
                    sn(is1) = sn(is1) + 1
                    goto 7100
                endif
            enddo
            ns = ns + 1
            srsn(ns) = a_rsn(1,isp1,ia1)
            !write(*,*) 'IAM\PSTA: ns = ', ns
7100 continue      
        endif
    enddo ! ia1
    enddo ! isp1
else
    ns = -9
endif

!
! population statisics
!
naa = 0
naf = 0
SRT = 0.

nda = 0.
cna = 0.
do nt1 = 1, 4
    do nt2 = 1, 4
        nma(nt1,nt2) = 0.
        nra(nt1,nt2) = 0.
    enddo
enddo
nmTa = 0.
nmXa = 0.
nmSa = 0.
nmXSa = 0.
nmNa = 0.
nrTa = 0.
nrXa = 0.
nrSa = 0.
nrXSa = 0.
nrNa = 0.
nmNmin = 999.
nmNmax = -999.
do nt1 = 1, 4
    fnta(nt1) = 0.
    fntta(nt1) = 0.
enddo
do iaa1 = 1, 20
    faaa(iaa1) = 0.
    do iaa2 = 1, 20
        naaxa(iaa1,iaa2) = 0.
    enddo
enddo
VmaxCa = 0.
VmaxNa = 0.
VmaxPa = 0.
VCa = 0.
VNa = 0.
VPa = 0.
q0DNACa = 0.
q0DNANa = 0.
q0DNAPa = 0.
q0aaCa = 0.
q0aaNa = 0.
qCa = 0.
qNa = 0.
qPa = 0.
CLimita = 0.
NLimita = 0.
PLimita = 0.
tga = 0.
kga = 0.
kgmin = 999.
kgmax = -999.
ntnka = 0.
ntnSka = 0.
ntnNka = 0.
ntnXka = 0.
ntnXkax = 0.

do isp1 = 1, dimnsp
naa_sp(isp1) = 0
do ia1 = 1, nat_sp(1,isp1)
    if (a_sa(1,isp1,ia1).ge.1) then
        ig1 = a_gid(1,isp1,ia1)
        naa_sp(isp1) = naa_sp(isp1) + 1
        naa = naa + 1
        SRT = SRT + a_sr(1,isp1,ia1)
        do nt1 = 1, 4
            fnta(nt1) = fnta(nt1) + a_ntnk(1,isp1,ia1,nt1) / dble(dimnnt) * a_sr(1,isp1,ia1)
        enddo
        do iaa1 = 1, 20
            faaa(iaa1) = faaa(iaa1) + dble(a_naa(1,isp1,ia1,iaa1)) / dble(g_naat(ig1)) * a_sr(1,isp1,ia1)
            do iaa2 = 1, 20
                naaxa(iaa1,iaa2) = naaxa(iaa1,iaa2) + dble(a_naax(1,isp1,ia1,iaa1,iaa2)) * a_sr(1,isp1,ia1)
            enddo
        enddo
        VmaxCa = VmaxCa + a_VmaxC(1,isp1,ia1) * a_sr(1,isp1,ia1)
        VmaxNa = VmaxNa + a_VmaxN(1,isp1,ia1) * a_sr(1,isp1,ia1)
        VmaxPa = VmaxPa + a_VmaxP(1,isp1,ia1) * a_sr(1,isp1,ia1)
        q0DNACa = q0DNACa + a_q0DNAC(1,isp1,ia1) * a_sr(1,isp1,ia1)
        q0DNANa = q0DNANa + a_q0DNAN(1,isp1,ia1) * a_sr(1,isp1,ia1)
        q0DNAPa = q0DNAPa + a_q0DNAP(1,isp1,ia1) * a_sr(1,isp1,ia1)
        q0aaCa = q0aaCa + a_q0aaC(1,isp1,ia1) * a_sr(1,isp1,ia1)
        q0aaNa = q0aaNa + a_q0aaN(1,isp1,ia1) * a_sr(1,isp1,ia1)
        if (oSlim.eq.0) then
            nda = nda + dble(a_nd(1,isp1,ia1)) * a_sr(1,isp1,ia1)
            do nt1 = 1, 4
                ntnka(nt1) = ntnka(nt1) + dble(a_ntnk(1,isp1,ia1,nt1)) * a_sr(1,isp1,ia1)
                ntnSka(nt1) = ntnSka(nt1) + dble(a_ntnSk(1,isp1,ia1,nt1)) * a_sr(1,isp1,ia1)
                ntnNka(nt1) = ntnNka(nt1) + dble(a_ntnNk(1,isp1,ia1,nt1)) * a_sr(1,isp1,ia1)
                ntnXka(nt1) = ntnXka(nt1) + dble(a_ntnXk(1,isp1,ia1,nt1)) * a_sr(1,isp1,ia1)
                ntnXkax(nt1) = ntnXkax(nt1) + dble(a_ntnXk(1,isp1,ia1,nt1))
            enddo
            cna = cna + dble(a_cn(1,isp1,ia1)) * a_sr(1,isp1,ia1)
            nmTa = nmTa + dble(a_nmT(1,isp1,ia1)) * a_sr(1,isp1,ia1)
            a_nmXt = a_nmT(1,isp1,ia1) - a_nmS(1,isp1,ia1) - a_nmN(1,isp1,ia1)
            a_nmXSt = a_nmXt + a_nmS(1,isp1,ia1)
            nmXa = nmXa + dble(a_nmXt) * a_sr(1,isp1,ia1)
            nmSa = nmSa + dble(a_nmS(1,isp1,ia1)) * a_sr(1,isp1,ia1)
            nmXSa = nmXSa + dble(a_nmXSt) * a_sr(1,isp1,ia1)
            nmNa = nmNa + dble(a_nmN(1,isp1,ia1)) * a_sr(1,isp1,ia1)
            nmNmin = min(a_nmN(1,isp1,ia1),nmNmin)
            nmNmax = max(a_nmN(1,isp1,ia1),nmNmax)
            nrta = nrta + dble(a_nrt(1,isp1,ia1)) * a_sr(1,isp1,ia1)
            nrSa = nrSa + dble(a_nrS(1,isp1,ia1)) * a_sr(1,isp1,ia1)
            nrNa = nrNa + dble(a_nrN(1,isp1,ia1)) * a_sr(1,isp1,ia1)
            do nt1 = 1, 4
                do nt2 = 1, 4
                    nma(nt1,nt2) = nma(nt1,nt2) + dble(a_nm(1,isp1,ia1,nt1,nt2)) * a_sr(1,isp1,ia1)
                    nra(nt1,nt2) = nra(nt1,nt2) + dble(a_nr(1,isp1,ia1,nt1,nt2)) * a_sr(1,isp1,ia1)
                    !if (nt1.ne.nt2) then
                    !    nmta = nmta + dble(a_nm(1,isp1,ia1,nt1,nt2)) * a_sr(1,isp1,ia1)
                    !    nrta = nrta + dble(a_nr(1,isp1,ia1,nt1,nt2)) * a_sr(1,isp1,ia1)
                    !endif
                enddo
            enddo
            do nt1 = 1, 4
                fntta(nt1) = fntta(nt1) + dble(a_nntt(1,isp1,ia1,nt1)) / g_nnttt(ig1) * a_sr(1,isp1,ia1)
            enddo
            VCa = VCa + a_VC(1,isp1,ia1) * a_sr(1,isp1,ia1)
            VNa = VNa + a_VN(1,isp1,ia1) * a_sr(1,isp1,ia1)
            VPa = VPa + a_VP(1,isp1,ia1) * a_sr(1,isp1,ia1)
            qCa = qCa + a_qC(1,isp1,ia1) * a_sr(1,isp1,ia1)
            qNa = qNa + a_qN(1,isp1,ia1) * a_sr(1,isp1,ia1)
            qPa = qPa + a_qP(1,isp1,ia1) * a_sr(1,isp1,ia1)
            CLimita = CLimita + a_CLimit(1,isp1,ia1) * a_sr(1,isp1,ia1)
            NLimita = NLimita + a_NLimit(1,isp1,ia1) * a_sr(1,isp1,ia1)
            PLimita = PLimita + a_PLimit(1,isp1,ia1) * a_sr(1,isp1,ia1)
            tga = tga + a_tg(1,isp1,ia1) * a_sr(1,isp1,ia1)
            
            if (t.eq.0.) then ! if t = 0, have to calculate kg
                if (oLN.eq.1) then
                    LNC = SC_sp(isp1) / (KmC + SC_sp(isp1))
                    LNN = SN_sp(isp1) / (KmN + SN_sp(isp1))
                    LNP = SP_sp(isp1) / (KmP + SP_sp(isp1))
                else
                    LNC = 1.
                    LNN = 1.
                    LNP = 1.
                endif
                call calc_V(isp1,ia1,LNC,LNN,LNP)
                call calc_q0(isp1,ia1,q0C,q0N,q0P)
                call calc_A(isp1,ia1,q0C,q0N,q0P,kg)
                if (kg.gt.0.) then
                    a_tg(1,isp1,ia1) = log(2.) / kg
                else
                    a_tg(1,isp1,ia1) = 1.e9
                endif
            endif
            if (a_tg(1,isp1,ia1).gt.0.) then
                kg = log(2.) / a_tg(1,isp1,ia1)
                kga = kga + kg * a_sr(1,isp1,ia1)
                kgmin = min(kg,kgmin)
                kgmax = max(kg,kgmax)
            endif
        endif
    endif ! sa
enddo ! ia1
naf = naf + naf_sp(isp1)
if (naa_sp(isp1).ne.(nat_sp(1,isp1)-naf_sp(isp1))) then
    write(*,*) 'IAM\ERR: naa_sp <> nat_sp-naf_sp.'
    write(*,*) 'IAM\ERR: naa_sp = ', naa_sp(isp1)
    write(*,*) 'IAM\ERR: nat_sp = ', nat_sp(1,isp1)
    write(*,*) 'IAM\ERR: naf_sp = ', naf_sp(isp1)
    goto 9100
endif
enddo ! isp1

if (isnan(SRT)) then
    write(*,*) 'IAM\ERR: SRT = nan.'
    goto 9100
endif

if (oHANS.eq.0) then
    if (naa.ne.int(SRT)) then
        write(*,*) 'IAM\ERR: naa <> SRT.'
        goto 9100
    endif
endif

!write(*,*) 'a2:', ntnXka(2)
!write(*,*) 'a3:', ntnXka(3)
!write(*,*) 'a4:', ntnXka(4)
!write(*,*) 'b1:', ntnXkax(1)
!write(*,*) 'b2:', ntnXkax(2)
!write(*,*) 'b3:', ntnXkax(3)
!write(*,*) 'b4:', ntnXkax(4)
!write(*,*) 'aSR:', SRT

nda = nda / SRT
do nt1 = 1, 4
    ntnka(nt1) = ntnka(nt1) / SRT
    ntnSka(nt1) = ntnSka(nt1) / SRT
    ntnNka(nt1) = ntnNka(nt1) / SRT
    ntnXka(nt1) = ntnXka(nt1) / SRT
enddo
cna = cna / SRT
do nt1 = 1, 4
    do nt2 = 1, 4
        nma(nt1,nt2) = nma(nt1,nt2) / SRT
        nra(nt1,nt2) = nra(nt1,nt2) / SRT
    enddo
enddo
nmTa = nmTa / SRT
nmXa = nmXa / SRT
nmSa = nmSa / SRT
nmXSa = nmXSa / SRT
nmNa = nmNa / SRT
nrta = nrta / SRT
nrSa = nrSa / SRT
nrNa = nrNa / SRT
!write(*,*) 'IAM\STA: DNA composition (fraction):'
do nt1 = 1, 4
    fnta(nt1) = fnta(nt1) / SRT
!    write(*,*) 'IAM\STA: ', nt1, fnta(nt1)
enddo
!write(*,*) 'IAM\STA: mRNA composition (fraction):'
do nt1 = 1, 4
    fntta(nt1) = fntta(nt1) / SRT
!    write(*,*) 'IAM\STA: ', nt1, fntta(nt1)
enddo
!write(*,*) 'IAM\STA: Amino acid composition (fraction):'
do iaa1 = 1, 20
    faaa(iaa1) = faaa(iaa1) / SRT
!    write(*,*) 'IAM\STA: ', iaa1, int_saa(iaa1), faaa(iaa1)
    do iaa2 = 1, 20
        naaxa(iaa1,iaa2) = naaxa(iaa1,iaa2) / SRT
    enddo
enddo
VCa = VCa / SRT
VNa = VNa / SRT
VPa = VPa / SRT
VmaxCa = VmaxCa / SRT
VmaxNa = VmaxNa / SRT
VmaxPa = VmaxPa / SRT
q0DNACa = q0DNACa / SRT
q0DNANa = q0DNANa / SRT
q0DNAPa = q0DNAPa / SRT
q0aaCa = q0aaCa / SRT
q0aaNa = q0aaNa / SRT
q0Ca = q0othC + q0DNACa + q0aaCa
q0Na = q0othN + q0DNANa + q0aaNa
q0Pa = q0othP + q0DNAPa 
qCa = qCa / SRT
qNa = qNa / SRT
qPa = qPa / SRT
CLimita = CLimita / SRT
NLimita = NLimita / SRT
PLimita = PLimita / SRT
tga = tga / SRT
kga = kga / SRT

! GC
ig1 = 1
GCT = (ntnka(3)+ntnka(4))/g_ntn(ig1)
GCX = (ntnXka(3)+ntnXka(4))/g_ntnX(ig1)
GCS = (ntnSka(3)+ntnSka(4))/g_ntnS(ig1)
GCXS = (ntnXka(3)+ntnXka(4)+ntnSka(3)+ntnSka(4))/(g_ntnX(ig1)+g_ntnS(ig1))
GCN = (ntnNka(3)+ntnNka(4))/g_ntnN(ig1)
write(*,*) 'IAM: GCT  = ', GCT
write(*,*) 'IAM: GCX  = ', GCX
write(*,*) 'IAM: GCS  = ', GCS
write(*,*) 'IAM: GCXS = ', GCXS
write(*,*) 'IAM: GCN  = ', GCN

! standard eviations
nmNstd = 0.
kgstd = 0.
do iaa1 = 1, 20
    faastd(iaa1) = 0.
enddo
if (oSlim.eq.0) then
    do isp1 = 1, dimnsp
        do ia1 = 1, nat_sp(1,isp1)
            if (a_sa(1,isp1,ia1).ge.1) then
                nmNstd = nmNstd + (a_nmN(1,isp1,ia1) - nmNa)**2 * a_sr(1,isp1,ia1)
                kgstd = kgstd + (kg - kga)**2 * a_sr(1,isp1,ia1)
                do iaa1 = 1, 20
                    faastd(iaa1) = faastd(iaa1) + (dble(a_naa(1,isp1,ia1,iaa1)) / dble(g_naat(ig1)) - faaa(iaa1))**2 * a_sr(1,isp1,ia1)
                enddo
            endif
        enddo
    enddo
    nmNstd = (nmNstd / SRT)**.5
    kgstd = (kgstd / SRT)**.5
    do iaa1 = 1, 20
        faastd(iaa1) = (faastd(iaa1) / SRT)**.5
    enddo
    if (ipSTAn.eq.1) then
        write(*,*) 'IAM\STA: Amino acid composition (fraction), mean, std:'
        do iaa1 = 1, 20
            write(*,*) 'IAM\STA: ', iaa1, int_saa(iaa1), faaa(iaa1), faastd(iaa1)
        enddo
    endif
endif

! correlations & anomalies
if (oSlim.eq.0) then
if (ipSTAn.eq.1) then
ig1 = 1

CT_nmT_nmT = 0.
CT_nmT_nmX = 0.
CT_nmT_nmS = 0.
CT_nmT_nmXS = 0.
CT_nmT_nmN = 0.
CT_nmX_nmT = 0.
CT_nmX_nmX = 0.
CT_nmX_nmS = 0.
CT_nmX_nmXS = 0.
CT_nmX_nmN = 0.
CT_nmS_nmT = 0.
CT_nmS_nmX = 0.
CT_nmS_nmS = 0.
CT_nmS_nmXS = 0.
CT_nmS_nmN = 0.
CT_nmXS_nmT = 0.
CT_nmXS_nmX = 0.
CT_nmXS_nmS = 0.
CT_nmXS_nmXS = 0.
CT_nmXS_nmN = 0.
CT_nmN_nmT = 0.
CT_nmN_nmX = 0.
CT_nmN_nmS = 0.
CT_nmN_nmXS = 0.
CT_nmN_nmN = 0.

CT_nmT = 0.
CT_nmX = 0.
CT_nmS = 0.
CT_nmXS = 0.
CT_nmN = 0.

CT_GCT_GCT = 0.
CT_GCT_GCX = 0.
CT_GCT_GCS = 0.
CT_GCT_GCXS = 0.
CT_GCT_GCN = 0.
CT_GCX_GCT = 0.
CT_GCX_GCX = 0.
CT_GCX_GCS = 0.
CT_GCX_GCXS = 0.
CT_GCX_GCN = 0.
CT_GCS_GCT = 0.
CT_GCS_GCX = 0.
CT_GCS_GCS = 0.
CT_GCS_GCXS = 0.
CT_GCS_GCN = 0.
CT_GCXS_GCT = 0.
CT_GCXS_GCX = 0.
CT_GCXS_GCS = 0.
CT_GCXS_GCXS = 0.
CT_GCXS_GCN = 0.
CT_GCN_GCT = 0.
CT_GCN_GCX = 0.
CT_GCN_GCS = 0.
CT_GCN_GCXS = 0.
CT_GCN_GCN = 0.

CT_GCT = 0.
CT_GCX = 0.
CT_GCS = 0.
CT_GCXS = 0.
CT_GCN = 0.

AT_nmT = 0.
AT_nmX = 0.
AT_nmS = 0.
AT_nmXS = 0.
AT_nmN = 0.

AT_GCT = 0.
AT_GCX = 0.
AT_GCS = 0.
AT_GCXS = 0.
AT_GCN = 0.

do isp1 = 1, dimnsp
do ia1 = 1, nat_sp(1,isp1)
    if (a_sa(1,isp1,ia1).ge.1) then

        a_nmXt = a_nmT(1,isp1,ia1) - a_nmS(1,isp1,ia1) - a_nmN(1,isp1,ia1)
        a_nmXSt = a_nmXt + a_nmS(1,isp1,ia1)
        
        a_GCTt = dble(a_ntnk(1,isp1,ia1,3)+a_ntnk(1,isp1,ia1,4))/g_ntn(ig1)
        a_GCXt = dble(a_ntnXk(1,isp1,ia1,3)+a_ntnXk(1,isp1,ia1,4))/g_ntnX(ig1)
        a_GCSt = dble(a_ntnSk(1,isp1,ia1,3)+a_ntnSk(1,isp1,ia1,4))/g_ntnS(ig1)
        a_GCXSt = dble(a_ntnXk(1,isp1,ia1,3)+a_ntnXk(1,isp1,ia1,4)+a_ntnSk(1,isp1,ia1,3)+a_ntnSk(1,isp1,ia1,4))/(g_ntnX(ig1)+g_ntnS(ig1))
        a_GCNt = dble(a_ntnNk(1,isp1,ia1,3)+a_ntnNk(1,isp1,ia1,4))/g_ntnN(ig1)

        CT_nmT_nmT = CT_nmT_nmT + (dble(a_nmT(1,isp1,ia1))-nmTa)*(dble(a_nmT(1,isp1,ia1))-nmTa) * a_sr(1,isp1,ia1)
        CT_nmT_nmX = CT_nmT_nmX + (dble(a_nmT(1,isp1,ia1))-nmTa)*(dble(a_nmXt)-nmXa) * a_sr(1,isp1,ia1)
        CT_nmT_nmS = CT_nmT_nmS + (dble(a_nmT(1,isp1,ia1))-nmTa)*(dble(a_nmS(1,isp1,ia1))-nmSa) * a_sr(1,isp1,ia1)
        CT_nmT_nmXS = CT_nmT_nmXS + (dble(a_nmT(1,isp1,ia1))-nmTa)*(dble(a_nmXSt)-nmXSa) * a_sr(1,isp1,ia1)
        CT_nmT_nmN = CT_nmT_nmN + (dble(a_nmT(1,isp1,ia1))-nmTa)*(dble(a_nmN(1,isp1,ia1))-nmNa) * a_sr(1,isp1,ia1)
        CT_nmX_nmT = CT_nmX_nmT + (dble(a_nmXt)-nmXa)*(dble(a_nmT(1,isp1,ia1))-nmTa) * a_sr(1,isp1,ia1)
        CT_nmX_nmX = CT_nmX_nmX + (dble(a_nmXt)-nmXa)*(dble(a_nmXt)-nmXa) * a_sr(1,isp1,ia1)
        CT_nmX_nmS = CT_nmX_nmS + (dble(a_nmXt)-nmXa)*(dble(a_nmS(1,isp1,ia1))-nmSa) * a_sr(1,isp1,ia1)
        CT_nmX_nmXS = CT_nmX_nmXS + (dble(a_nmXt)-nmXa)*(dble(a_nmXSt)-nmXSa) * a_sr(1,isp1,ia1)
        CT_nmX_nmN = CT_nmX_nmN + (dble(a_nmXt)-nmXa)*(dble(a_nmN(1,isp1,ia1))-nmNa) * a_sr(1,isp1,ia1)
        CT_nmS_nmT = CT_nmS_nmT + (dble(a_nmS(1,isp1,ia1))-nmSa)*(dble(a_nmT(1,isp1,ia1))-nmTa) * a_sr(1,isp1,ia1)
        CT_nmS_nmX = CT_nmS_nmX + (dble(a_nmS(1,isp1,ia1))-nmSa)*(dble(a_nmXt)-nmXa) * a_sr(1,isp1,ia1)
        CT_nmS_nmS = CT_nmS_nmS + (dble(a_nmS(1,isp1,ia1))-nmSa)*(dble(a_nmS(1,isp1,ia1))-nmSa) * a_sr(1,isp1,ia1)
        CT_nmS_nmXS = CT_nmS_nmXS + (dble(a_nmS(1,isp1,ia1))-nmSa)*(dble(a_nmXSt)-nmXSa) * a_sr(1,isp1,ia1)
        CT_nmS_nmN = CT_nmS_nmN + (dble(a_nmS(1,isp1,ia1))-nmSa)*(dble(a_nmN(1,isp1,ia1))-nmNa) * a_sr(1,isp1,ia1)
        CT_nmXS_nmT = CT_nmXS_nmT + (dble(a_nmXSt)-nmXSa)*(dble(a_nmT(1,isp1,ia1))-nmTa) * a_sr(1,isp1,ia1)
        CT_nmXS_nmX = CT_nmXS_nmX + (dble(a_nmXSt)-nmXSa)*(dble(a_nmXt)-nmXa) * a_sr(1,isp1,ia1)
        CT_nmXS_nmS = CT_nmXS_nmS + (dble(a_nmXSt)-nmXSa)*(dble(a_nmS(1,isp1,ia1))-nmSa) * a_sr(1,isp1,ia1)
        CT_nmXS_nmXS = CT_nmXS_nmXS + (dble(a_nmXSt)-nmXSa)*(dble(a_nmXSt)-nmXSa) * a_sr(1,isp1,ia1)
        CT_nmXS_nmN = CT_nmXS_nmN + (dble(a_nmXSt)-nmXSa)*(dble(a_nmN(1,isp1,ia1))-nmNa) * a_sr(1,isp1,ia1)
        CT_nmN_nmT = CT_nmN_nmT + (dble(a_nmN(1,isp1,ia1))-nmNa)*(dble(a_nmT(1,isp1,ia1))-nmTa) * a_sr(1,isp1,ia1)
        CT_nmN_nmX = CT_nmN_nmX + (dble(a_nmN(1,isp1,ia1))-nmNa)*(dble(a_nmXt)-nmXa) * a_sr(1,isp1,ia1)
        CT_nmN_nmS = CT_nmN_nmS + (dble(a_nmN(1,isp1,ia1))-nmNa)*(dble(a_nmS(1,isp1,ia1))-nmSa) * a_sr(1,isp1,ia1)
        CT_nmN_nmXS = CT_nmN_nmXS + (dble(a_nmN(1,isp1,ia1))-nmNa)*(dble(a_nmXSt)-nmXSa) * a_sr(1,isp1,ia1)
        CT_nmN_nmN = CT_nmN_nmN + (dble(a_nmN(1,isp1,ia1))-nmNa)*(dble(a_nmN(1,isp1,ia1))-nmNa) * a_sr(1,isp1,ia1)

        CT_nmT = CT_nmT + (dble(a_nmT(1,isp1,ia1))-nmTa)**2 * a_sr(1,isp1,ia1)
        CT_nmX = CT_nmX + (dble(a_nmXt)-nmXa)**2 * a_sr(1,isp1,ia1)
        CT_nmS = CT_nmS + (dble(a_nmS(1,isp1,ia1))-nmSa)**2 * a_sr(1,isp1,ia1)
        CT_nmXS = CT_nmXS + (dble(a_nmXSt)-nmXa-nmSa)**2 * a_sr(1,isp1,ia1)
        CT_nmN = CT_nmN + (dble(a_nmN(1,isp1,ia1))-nmNa)**2 * a_sr(1,isp1,ia1)

        CT_GCT_GCT = CT_GCT_GCT + (a_GCTt-GCT)*(a_GCTt-GCT) * a_sr(1,isp1,ia1)
        CT_GCT_GCX = CT_GCT_GCX + (a_GCTt-GCT)*(a_GCXt-GCX) * a_sr(1,isp1,ia1)
        CT_GCT_GCS = CT_GCT_GCS + (a_GCTt-GCT)*(a_GCSt-GCS) * a_sr(1,isp1,ia1)
        CT_GCT_GCXS = CT_GCT_GCXS + (a_GCTt-GCT)*(a_GCXSt-GCXS) * a_sr(1,isp1,ia1)
        CT_GCT_GCN = CT_GCT_GCN + (a_GCTt-GCT)*(a_GCNt-GCN) * a_sr(1,isp1,ia1)
        CT_GCX_GCT = CT_GCX_GCT + (a_GCXt-GCX)*(a_GCTt-GCT) * a_sr(1,isp1,ia1)
        CT_GCX_GCX = CT_GCX_GCX + (a_GCXt-GCX)*(a_GCXt-GCX) * a_sr(1,isp1,ia1)
        CT_GCX_GCS = CT_GCX_GCS + (a_GCXt-GCX)*(a_GCSt-GCS) * a_sr(1,isp1,ia1)
        CT_GCX_GCXS = CT_GCX_GCXS + (a_GCXt-GCX)*(a_GCXSt-GCXS) * a_sr(1,isp1,ia1)
        CT_GCX_GCN = CT_GCX_GCN + (a_GCXt-GCX)*(a_GCNt-GCN) * a_sr(1,isp1,ia1)
        CT_GCS_GCT = CT_GCS_GCT + (a_GCSt-GCS)*(a_GCTt-GCT) * a_sr(1,isp1,ia1)
        CT_GCS_GCX = CT_GCS_GCX + (a_GCSt-GCS)*(a_GCXt-GCX) * a_sr(1,isp1,ia1)
        CT_GCS_GCS = CT_GCS_GCS + (a_GCSt-GCS)*(a_GCSt-GCS) * a_sr(1,isp1,ia1)
        CT_GCS_GCXS = CT_GCS_GCXS + (a_GCSt-GCS)*(a_GCXSt-GCXS) * a_sr(1,isp1,ia1)
        CT_GCS_GCN = CT_GCS_GCN + (a_GCSt-GCS)*(a_GCNt-GCN) * a_sr(1,isp1,ia1)
        CT_GCXS_GCT = CT_GCXS_GCT + (a_GCXSt-GCXS)*(a_GCTt-GCT) * a_sr(1,isp1,ia1)
        CT_GCXS_GCX = CT_GCXS_GCX + (a_GCXSt-GCXS)*(a_GCXt-GCX) * a_sr(1,isp1,ia1)
        CT_GCXS_GCS = CT_GCXS_GCS + (a_GCXSt-GCXS)*(a_GCSt-GCS) * a_sr(1,isp1,ia1)
        CT_GCXS_GCXS = CT_GCXS_GCXS + (a_GCXSt-GCXS)*(a_GCXSt-GCXS) * a_sr(1,isp1,ia1)
        CT_GCXS_GCN = CT_GCXS_GCN + (a_GCXSt-GCXS)*(a_GCNt-GCN) * a_sr(1,isp1,ia1)
        CT_GCN_GCT = CT_GCN_GCT + (a_GCNt-GCN)*(a_GCTt-GCT) * a_sr(1,isp1,ia1)
        CT_GCN_GCX = CT_GCN_GCX + (a_GCNt-GCN)*(a_GCXt-GCX) * a_sr(1,isp1,ia1)
        CT_GCN_GCS = CT_GCN_GCS + (a_GCNt-GCN)*(a_GCSt-GCS) * a_sr(1,isp1,ia1)
        CT_GCN_GCXS = CT_GCN_GCXS + (a_GCNt-GCN)*(a_GCXSt-GCXS) * a_sr(1,isp1,ia1)
        CT_GCN_GCN = CT_GCN_GCN + (a_GCNt-GCN)*(a_GCNt-GCN) * a_sr(1,isp1,ia1)

        CT_GCT = CT_GCT + (a_GCTt-GCT)**2 * a_sr(1,isp1,ia1)
        CT_GCX = CT_GCX + (a_GCXt-GCX)**2 * a_sr(1,isp1,ia1)
        CT_GCS = CT_GCS + (a_GCSt-GCS)**2 * a_sr(1,isp1,ia1)
        CT_GCXS = CT_GCXS + (a_GCXSt-GCXS)**2 * a_sr(1,isp1,ia1)
        CT_GCN = CT_GCN + (a_GCNt-GCN)**2 * a_sr(1,isp1,ia1)

        AT_nmT = AT_nmT + (dble(a_nmT(1,isp1,ia1))-0.) * a_sr(1,isp1,ia1)
        AT_nmX = AT_nmX + (dble(a_nmXt)-0.) * a_sr(1,isp1,ia1)
        AT_nmS = AT_nmS + (dble(a_nmS(1,isp1,ia1))-0.) * a_sr(1,isp1,ia1)
        AT_nmXS = AT_nmXS + (dble(a_nmXSt)-0.-0.) * a_sr(1,isp1,ia1)
        AT_nmN = AT_nmN + (dble(a_nmN(1,isp1,ia1))-0.) * a_sr(1,isp1,ia1)

        AT_GCT = AT_GCT + (a_GCTt-(g_ntnk(ig1,3)+g_ntnk(ig1,4)) / g_ntn(ig1)) * a_sr(1,isp1,ia1)
        AT_GCX = AT_GCX + (a_GCXt-(g_ntnXk(ig1,3)+g_ntnXk(ig1,4)) / g_ntnX(ig1)) * a_sr(1,isp1,ia1)
        AT_GCS = AT_GCS + (a_GCSt-(g_ntnSk(ig1,3)+g_ntnSk(ig1,4)) / g_ntnS(ig1)) * a_sr(1,isp1,ia1)
        AT_GCXS = AT_GCXS + (a_GCXSt-((g_ntnXk(ig1,3)+g_ntnXk(ig1,4))+(g_ntnSk(ig1,3)+g_ntnSk(ig1,4))) / (g_ntnX(ig1)+g_ntnS(ig1))) * a_sr(1,isp1,ia1)
        AT_GCN = AT_GCN + (a_GCNt-(g_ntnNk(ig1,3)+g_ntnNk(ig1,4)) / g_ntnN(ig1)) * a_sr(1,isp1,ia1)

        
    endif
enddo
enddo
write(*,*) 'IAM: Correlation coefficients.'
write(*,*) ''
write(*,*) 'IAM: nmT vs. nmT  = ', CT_nmT_nmT / (CT_nmT*CT_nmT)**.5
write(*,*) 'IAM: nmT vs. nmX  = ', CT_nmT_nmX / (CT_nmT*CT_nmX)**.5
write(*,*) 'IAM: nmT vs. nmS  = ', CT_nmT_nmS / (CT_nmT*CT_nmS)**.5
write(*,*) 'IAM: nmT vs. nmSX = ', CT_nmT_nmXS / (CT_nmT*CT_nmXS)**.5
write(*,*) 'IAM: nmT vs. nmN  = ', CT_nmT_nmN / (CT_nmT*CT_nmN)**.5
write(*,*) ''
write(*,*) 'IAM: nmX vs. nmT  = ', CT_nmX_nmT / (CT_nmX*CT_nmT)**.5
write(*,*) 'IAM: nmX vs. nmX  = ', CT_nmX_nmX / (CT_nmX*CT_nmX)**.5
write(*,*) 'IAM: nmX vs. nmS  = ', CT_nmX_nmS / (CT_nmX*CT_nmS)**.5
write(*,*) 'IAM: nmX vs. nmSX = ', CT_nmX_nmXS / (CT_nmX*CT_nmXS)**.5
write(*,*) 'IAM: nmX vs. nmN  = ', CT_nmX_nmN / (CT_nmX*CT_nmN)**.5
write(*,*) ''
write(*,*) 'IAM: nmS vs. nmT  = ', CT_nmS_nmT / (CT_nmS*CT_nmT)**.5
write(*,*) 'IAM: nmS vs. nmX  = ', CT_nmS_nmX / (CT_nmS*CT_nmX)**.5
write(*,*) 'IAM: nmS vs. nmS  = ', CT_nmS_nmS / (CT_nmS*CT_nmS)**.5
write(*,*) 'IAM: nmS vs. nmSX = ', CT_nmS_nmXS / (CT_nmS*CT_nmXS)**.5
write(*,*) 'IAM: nmS vs. nmN  = ', CT_nmS_nmN / (CT_nmS*CT_nmN)**.5
write(*,*) ''
write(*,*) 'IAM: nmXS vs. nmT  = ', CT_nmXS_nmT / (CT_nmXS*CT_nmT)**.5
write(*,*) 'IAM: nmXS vs. nmX  = ', CT_nmXS_nmX / (CT_nmXS*CT_nmX)**.5
write(*,*) 'IAM: nmXS vs. nmS  = ', CT_nmXS_nmS / (CT_nmXS*CT_nmS)**.5
write(*,*) 'IAM: nmXS vs. nmSX = ', CT_nmXS_nmXS / (CT_nmXS*CT_nmXS)**.5
write(*,*) 'IAM: nmXS vs. nmN  = ', CT_nmXS_nmN / (CT_nmXS*CT_nmN)**.5
write(*,*) ''
write(*,*) 'IAM: nmN vs. nmT  = ', CT_nmN_nmT / (CT_nmN*CT_nmT)**.5
write(*,*) 'IAM: nmN vs. nmX  = ', CT_nmN_nmX / (CT_nmN*CT_nmX)**.5
write(*,*) 'IAM: nmN vs. nmS  = ', CT_nmN_nmS / (CT_nmN*CT_nmS)**.5
write(*,*) 'IAM: nmN vs. nmSX = ', CT_nmN_nmXS / (CT_nmN*CT_nmXS)**.5
write(*,*) 'IAM: nmN vs. nmN  = ', CT_nmN_nmN / (CT_nmN*CT_nmN)**.5
write(*,*) ''
write(*,*) 'IAM: GCT vs. GCT  = ', CT_GCT_GCT / (CT_GCT*CT_GCT)**.5
write(*,*) 'IAM: GCT vs. GCX  = ', CT_GCT_GCX / (CT_GCT*CT_GCX)**.5
write(*,*) 'IAM: GCT vs. GCS  = ', CT_GCT_GCS / (CT_GCT*CT_GCS)**.5
write(*,*) 'IAM: GCT vs. GCXS = ', CT_GCT_GCXS / (CT_GCT*CT_GCXS)**.5
write(*,*) 'IAM: GCT vs. GCN  = ', CT_GCT_GCN / (CT_GCT*CT_GCN)**.5
write(*,*) ''
write(*,*) 'IAM: GCX vs. GCT  = ', CT_GCX_GCT / (CT_GCX*CT_GCT)**.5
write(*,*) 'IAM: GCX vs. GCX  = ', CT_GCX_GCX / (CT_GCX*CT_GCX)**.5
write(*,*) 'IAM: GCX vs. GCS  = ', CT_GCX_GCS / (CT_GCX*CT_GCS)**.5
write(*,*) 'IAM: GCX vs. GCXS = ', CT_GCX_GCXS / (CT_GCX*CT_GCXS)**.5
write(*,*) 'IAM: GCX vs. GCN  = ', CT_GCX_GCN / (CT_GCX*CT_GCN)**.5
write(*,*) ''
write(*,*) 'IAM: GCS vs. GCT  = ', CT_GCS_GCT / (CT_GCS*CT_GCT)**.5
write(*,*) 'IAM: GCS vs. GCX  = ', CT_GCS_GCX / (CT_GCS*CT_GCX)**.5
write(*,*) 'IAM: GCS vs. GCS  = ', CT_GCS_GCS / (CT_GCS*CT_GCS)**.5
write(*,*) 'IAM: GCS vs. GCXS = ', CT_GCS_GCXS / (CT_GCS*CT_GCXS)**.5
write(*,*) 'IAM: GCS vs. GCN  = ', CT_GCS_GCN / (CT_GCS*CT_GCN)**.5
write(*,*) ''
write(*,*) 'IAM: GCXS vs. GCT  = ', CT_GCXS_GCT / (CT_GCXS*CT_GCT)**.5
write(*,*) 'IAM: GCXS vs. GCX  = ', CT_GCXS_GCX / (CT_GCXS*CT_GCX)**.5
write(*,*) 'IAM: GCXS vs. GCS  = ', CT_GCXS_GCS / (CT_GCXS*CT_GCS)**.5
write(*,*) 'IAM: GCXS vs. GCXS = ', CT_GCXS_GCXS / (CT_GCXS*CT_GCXS)**.5
write(*,*) 'IAM: GCXS vs. GCN  = ', CT_GCXS_GCN / (CT_GCXS*CT_GCN)**.5
write(*,*) ''
write(*,*) 'IAM: GCN vs. GCT  = ', CT_GCN_GCT / (CT_GCN*CT_GCT)**.5
write(*,*) 'IAM: GCN vs. GCX  = ', CT_GCN_GCX / (CT_GCN*CT_GCX)**.5
write(*,*) 'IAM: GCN vs. GCS  = ', CT_GCN_GCS / (CT_GCN*CT_GCS)**.5
write(*,*) 'IAM: GCN vs. GCXS = ', CT_GCN_GCXS / (CT_GCN*CT_GCXS)**.5
write(*,*) 'IAM: GCN vs. GCN  = ', CT_GCN_GCN / (CT_GCN*CT_GCN)**.5
write(*,*) ''

write(*,*) 'IAM: Anomalies.'
write(*,*) ''
write(*,*) 'IAM: nmT  = ', AT_nmT / SRT
write(*,*) 'IAM: nmX  = ', AT_nmX / SRT
write(*,*) 'IAM: nmS  = ', AT_nmS / SRT
write(*,*) 'IAM: nmXS = ', AT_nmXS / SRT
write(*,*) 'IAM: nmN  = ', AT_nmN / SRT
write(*,*) ''
write(*,*) 'IAM: GCT  = ', AT_GCT / SRT
write(*,*) 'IAM: GCX  = ', AT_GCX / SRT
write(*,*) 'IAM: GCS  = ', AT_GCS / SRT
write(*,*) 'IAM: GCXS = ', AT_GCXS / SRT
write(*,*) 'IAM: GCN  = ', AT_GCN / SRT

endif
endif

if (rn.gt.0.) then
    rdelta = rdeltat / rn
    rnu = rnut / rn
else
    rdelta = -9.
    rnu = -9.
endif    

if (fpn2.gt.0.) then
    fp = fpt / fpn2
else
    fp = -9.
endif    


!   
! --- change statistics ---
!
write(*,*) 'IAM:'
write(*,*) 'IAM: Change statistics:'
ig1 = 1
t_now = t
GC_now = fnta(3)+fnta(4)
GCT_now = GCT
GCS_now = GCS
GCN_now = GCN
GCX_now = GCX
GCXS_now = GCXS
write(*,*) 'IAM: nmSa    = ', nmSa
write(*,*) 'IAM: nrSa    = ', nrSa
write(*,*) 'IAM: nmNa    = ', nmNa
write(*,*) 'IAM: nrNa    = ', nrNa
write(*,*) 'IAM: nmta    = ', nmta
write(*,*) 'IAM: nrta    = ', nrta
write(*,*) 'IAM: g_ncSf  = ', g_ncSf(ig1)
write(*,*) 'IAM: g_ncSff = ', g_ncSff(ig1)
write(*,*) 'IAM: g_ncNf  = ', g_ncNf(ig1)
write(*,*) 'IAM: g_ncNff = ', g_ncNff(ig1)
write(*,*) 'IAM: g_ncmf  = ', g_ncmf(ig1)
write(*,*) 'IAM: g_ncmff = ', g_ncmff(ig1)
write(*,*) 'IAM: g_ncrf  = ', g_ncrf(ig1)
write(*,*) 'IAM: g_ncrff = ', g_ncrff(ig1)
write(*,*) 'IAM: g_ntnS  = ', g_ntnS(ig1)
write(*,*) 'IAM: g_ntnN  = ', g_ntnN(ig1)
write(*,*) 'IAM: g_ntn   = ', g_ntn(ig1)
piS_now = (nmSa+nrSa-dble(g_ncSf(ig1)+g_ncSff(ig1)))/g_ntnS(ig1)
piN_now = (nmNa+nrNa-dble(g_ncNf(ig1)+g_ncNff(ig1)))/g_ntnN(ig1)
piM_now = (nmta-dble(g_ncmf(ig1)+g_ncmff(ig1)))/g_ntn(ig1)
piR_now = (nrta-dble(g_ncrf(ig1)+g_ncrff(ig1)))/g_ntn(ig1)
piT_now = piM_now+piR_now
dNdS_now = -9.
if ((nmNt.gt.0.).and.(nmSa.gt.0.).and.(nmSt.gt.0.)) then
    dNdS_now = (nmNa/g_ntnN(ig1))/(nmSa/g_ntnS(ig1))
endif
rm_now = -9.
if (nmta.gt.0.) then
    rm_now = nrta/nmta
endif
q0DNAC_now = q0DNACa
q0DNAN_now = q0DNANa
q0DNAP_now = q0DNAPa
q0aaC_now = q0aaCa
q0aaN_now = q0aaNa
q0C_now = q0Ca
q0N_now = q0Na
q0P_now = q0Pa
kg_now = kga
write(*,*) 'IAM: GC now       = ', GC_now
write(*,*) 'IAM: GCT now      = ', GCT_now
write(*,*) 'IAM: GCS now      = ', GCS_now
write(*,*) 'IAM: GCN now      = ', GCN_now
write(*,*) 'IAM: GCX now      = ', GCX_now
write(*,*) 'IAM: GCXS now     = ', GCXS_now

write(*,*) 'IAM: piS(low) now = ', piS_now
if (piS_now.lt.piS_last) then
    if (((piS_last-piS_now)/piS_last*100.).gt.1.) then
        write(*,*) 'IAM: piS decrease %, t(y) = ', int((piS_last-piS_now)/piS_last*100.), int(t/365.25)
    endif
endif
write(*,*) 'IAM: piN(low) now = ', piN_now
write(*,*) 'IAM: piM(low) now = ', piM_now
write(*,*) 'IAM: piR(low) now = ', piR_now
write(*,*) 'IAM: piT(low) now = ', piT_now
write(*,*) 'IAM: dNdS now     = ', dNdS_now
write(*,*) 'IAM: rm now       = ', rm_now
write(*,*) 'IAM: q0DNAC now   = ', q0DNAC_now
write(*,*) 'IAM: q0DNAN now   = ', q0DNAN_now
write(*,*) 'IAM: q0DNAP now   = ', q0DNAP_now
write(*,*) 'IAM: q0aaC now    = ', q0aaC_now
write(*,*) 'IAM: q0aaN now    = ', q0aaN_now
write(*,*) 'IAM: q0C now      = ', q0C_now
write(*,*) 'IAM: q0N now      = ', q0N_now
write(*,*) 'IAM: q0P now      = ', q0P_now
write(*,*) 'IAM: kg now       = ', kg_now

t_last = t_now
GC_last = GC_now
GCT_last = GCT_now
GCS_last = GCS_now
GCN_last = GCN_now
GCX_last = GCX_now
GCXS_last = GCXS_now
piS_last = piS_now
piN_last = piN_now
piM_last = piM_now
piR_last = piR_now
piT_last = piT_now
dNdS_last = dNdS_now
rm_last = rm_now
q0DNAC_last = q0DNAC_now
q0DNAN_last = q0DNAN_now
q0DNAP_last = q0DNAP_now
q0aaC_last = q0aaC_now
q0aaN_last = q0aaN_now
q0C_last = q0C_now
q0N_last = q0N_now
q0P_last = q0P_now
kg_last = kg_now
if (t_first.lt.0.) then
    t_first = t_now
    GC_first = GC_now
    GCT_first = GCT_now
    GCS_first = GCS_now
    GCN_first = GCN_now
    GCX_first = GCX_now
    GCXS_first = GCXS_now
    piS_first = piS_now
    piN_first = piN_now
    piM_first = piM_now
    piR_first = piR_now
    piT_first = piT_now
    dNdS_first = dNdS_now
    rm_first = rm_now
    q0DNAC_first = q0DNAC_now
    q0DNAN_first = q0DNAN_now
    q0DNAP_first = q0DNAP_now
    q0aaC_first = q0aaC_now
    q0aaN_first = q0aaN_now
    q0C_first = q0C_now
    q0N_first = q0N_now
    q0P_first = q0P_now
    kg_first = kg_now
endif

if ((ipSTAf.eq.1).and.(opSTAf.eq.0)) then
    write(*,*) 'IAM\PSTA: Final, skipping STA print.'
    goto 9100
endif

!
! --- write ---
!
isp1 = 1
ig1 = 1

write(241,241) t, dble(nat), dble(naa), dble(naf), dble(ns), SRT, LP, K, nadt, namt, nmSt, nmNt, nart, nrSt, nrNt, rdelta, rnu, &
    nawt, SC_sp(isp1), SN_sp(isp1), SP_sp(isp1), Q, fp, -9., th, cp_delta, g_ntnS(ig1), g_ntnN(ig1), -9., &
    nda, &
    nma(1,2), nma(1,3), nma(1,4), nma(2,1), nma(2,3), nma(2,4), nma(3,1), nma(3,2), nma(3,4), nma(4,1), nma(4,2), nma(4,3), &
    nra(1,2), nra(1,3), nra(1,4), nra(2,1), nra(2,3), nra(2,4), nra(3,1), nra(3,2), nra(3,4), nra(4,1), nra(4,2), nra(4,3), &
    nmSa, nmNa, nrSa, nrNa, dble(g_ncf(1)+g_ncff(1)), &
    dble(g_ncSf(1)+g_ncSff(1)), dble(g_ncNf(1)+g_ncNff(1)), dble(g_ncmf(1)+g_ncmff(1)), dble(g_ncrf(1)+g_ncrff(1)), nmNstd, &
    cna, nrta, &
    (fnta(nt1), nt1 = 1, 4), (fntta(nt1), nt1 = 1, 4), (faaa(iaa1), iaa1 = 1, 20), &
    VmaxCa, VmaxNa, VmaxPa, VCa, VNa, VPa, &
     q0DNACa, q0DNANa, q0DNAPa, q0aaCa, q0aaNa, q0Ca, q0Na, q0Pa, qCa, qNa, qPa, CLimita, NLimita, PLimita, &
      tga, kga, kgmin, kgmax, kgstd, &
    GCT, GCS, GCN, GCX, GCXS

241 format(124e24.16)

    do iaa1 = 1, 20
        write(242,242) t, dble(iaa1), (naaxa(iaa1,iaa2), iaa2 = 1, 20)
    enddo

242 format(22e24.16)
    
tpnSTA = tpnSTA + dtpSTA

if (ipSTAf.eq.1) then
    goto 9100
endif

endif

!    
! ------------------
! --- output GLO ---
! ------------------
!
if ((opGLO.gt.0).and.(t.ge.tpnGLO)) then

ipGLO = ipGLO + 1
write(*,*) 'IAM\GLO: Writing GLO genomes. ipGLO = ', ipGLO

if (cTime.gt.0) then
    open(221,file='IAM_O_DNA_GLO_'//sap(ipGLO)//cString//'.txt')
    open(222,file='IAM_O_VAR_GLO_'//sap(ipGLO)//cString//'.txt')
else
    open(221,file='IAM_O_DNA_GLO_'//sap(ipGLO)//'.txt')
    open(222,file='IAM_O_VAR_GLO_'//sap(ipGLO)//'.txt')
endif

if ((oHANS.ge.1).and.(opGLO.eq.1)) then
    SRT = 0.
    do isp1 = 1, dimnsp
        do ia1 = 1, nat_sp(1,isp1)
            if (a_sa(1,isp1,ia1).ge.1) then
                SRT = SRT + a_sr(1,isp1,ia1)
            endif
        enddo
    enddo
endif

isp1 = 1
ia2 = 0
do ia1 = 1, npGLO
7610 continue
    if (opGLO.eq.1) then
        if (oHANS.ge.1) then
            r = r4_uni(jdum)
            SRCDF = 0.
            do isp1 = 1, dimnsp
            do ia2 = 1, nat_sp(1,isp1)
            if (a_sa(1,isp1,ia2).ge.1) then
                SRCDF = SRCDF + a_sr(1,isp1,ia2) / SRT
                if (SRCDF.ge.r) then
                    goto 7620
                endif
            endif
            enddo
            enddo
            write(*,*) 'IAM\GLO: ERR on agent selection.'
            goto 9100
7620 continue
        else
            r = r4_uni(jdum)
            isp1 = 1 + r * dimnsp
            r = r4_uni(jdum)
            ia2 = 1 + r * nat_sp(1,isp1)
        endif
    else
        ia2 = ia2 + 1
        if (ia2.gt.nat_sp(1,isp1)) then
            isp1 = isp1 + 1
            ia2 = 1
        endif
        if (isp1.gt.dimnsp) then
            write(*,*) 'IAM\GLO: Warning: Not enough active agents.'
            goto 7630
        endif
    endif
    ia3 = ia2
    !write(*,*) 'IAM\GLO: ia1, ia2, ia3 = ', ia1, ia2, ia3
    if (a_sa(1,isp1,ia2).eq.0) then
        !write(*,*) 'IAM\GLO: Inactive agent, trying again...'
        goto 7610
    endif
    write(222,206) t, ia1, ia3, a_sa(1,isp1,ia3), a_sn(1,isp1,ia3), a_rsn(1,isp1,ia3), a_nd(1,isp1,ia3), &
        a_nm(1,isp1,ia3,1,2), a_nm(1,isp1,ia3,1,3), a_nm(1,isp1,ia3,1,4), a_nm(1,isp1,ia3,2,1), a_nm(1,isp1,ia3,2,3), a_nm(1,isp1,ia3,2,4), a_nm(1,isp1,ia3,3,1), a_nm(1,isp1,ia3,3,2), a_nm(1,isp1,ia3,3,4), a_nm(1,isp1,ia3,4,1), a_nm(1,isp1,ia3,4,2), a_nm(1,isp1,ia3,4,3), &
        -9, a_tb(1,isp1,ia3), -9, a_gid(1,isp1,ia3), a_cn(1,isp1,ia3)
    do int1 = 1, dimnnt
        g_ntu1(1,int1) = g_nt(a_gid(1,isp1,ia3),int1)
    enddo
    do ic1 = 1, a_cn(1,isp1,ia3)
        g_ntu1(1,a_ck(ic1,isp1,ia3,1)) = a_cnt(ic1,isp1,ia3,1)
    enddo
    !write(*,*) '>IAMGLO'//sap(ipGLO)//sap(ia1)//' IAM GENOME GLO ip: '//sap(ipGLO)//' ia: '//sap(ia1)//'.'
    write(221,205) '>IAMGLO'//sap(ipGLO)//sap(ia1)//' IAM GENOME GLO ip: '//sap(ipGLO)//' ia: '//sap(ia1)//'.'
    write(221,201) (int_snt(g_ntu1(1,int1)),int1=1,dimnnt)
enddo ! ia1

7630 continue
     
close(221)
close(222)

tpnGLO = tpnGLO + dtpGLO

endif  

!
! ------------------
! --- output CLK ---
! ------------------
!
if ((opCLK.eq.1).and.(ipCLK.eq.0)) then

ipCLK = 1
write(*,*) 'IAM\G: Writing CLK genome. ipCLK = ', ipCLK

if (cTime.gt.0) then
    open(231,file='IAM_O_DNA_CLK'//cString//'.txt')
    open(232,file='IAM_O_VAR_CLK'//cString//'.txt')
else
    open(231,file='IAM_O_DNA_CLK.txt')
    open(232,file='IAM_O_VAR_CLK.txt')
endif

ia1 = 1
r = r4_uni(jdum)
isp1 = 1 + r * dimnsp
r = r4_uni(jdum)
ia2 = 1 + r * nat_sp(1,isp1)
ia3 = ia2
!write(*,*) 'IAM\CLK: ia1, ia2, ia3 = ', ia1, ia2, ia3
write(232,206) t, ia1, ia3, a_sa(1,isp1,ia3), a_sn(1,isp1,ia3), a_rsn(1,isp1,ia3), a_nd(1,isp1,ia3), &
    a_nm(1,isp1,ia3,1,2), a_nm(1,isp1,ia3,1,3), a_nm(1,isp1,ia3,1,4), a_nm(1,isp1,ia3,2,1), a_nm(1,isp1,ia3,2,3), a_nm(1,isp1,ia3,2,4), a_nm(1,isp1,ia3,3,1), a_nm(1,isp1,ia3,3,2), a_nm(1,isp1,ia3,3,4), a_nm(1,isp1,ia3,4,1), a_nm(1,isp1,ia3,4,2), a_nm(1,isp1,ia3,4,3), &
    -9, a_tb(1,isp1,ia3), -9, a_gid(1,isp1,ia3), a_cn(1,isp1,ia3)
do int1 = 1, dimnnt
    g_ntu1(1,int1) = g_nt(a_gid(1,isp1,ia3),int1)
enddo
!write(*,*) '>IAMCLK'//sap(ipCLK)//sap(ia1)//' IAM GENOME CLK ip: '//sap(ipCLK)//' ia: '//sap(ia1)//'.'
write(231,205) '>IAMCLK'//sap(ipCLK)//sap(ia1)//' IAM GENOME CLK ip: '//sap(ipCLK)//' ia: '//sap(ia1)//'.'
write(231,201) (int_snt(g_ntu1(1,int1)),int1=1,dimnnt)

close(231)
close(232)

endif  

!
! ------------------
! --- output POP ---
! ------------------
!
if (t.ge.tpnPOP) then
ipPOP = ipPOP + 1
write(*,*) 'IAM\PPOP: t, ipPOP = ', t, ipPOP

if (cTime.gt.0) then
    open(271,file='IAM_O_POP_'//sap(ipPOP)//cString//'.txt')
else
    open(271,file='IAM_O_POP_'//sap(ipPOP)//'.txt')
endif

do isp1 = 1, dimnsp
do ia1 = 1, nat_sp(1,isp1)
    if (a_sa(1,isp1,ia1).ge.1) then

        if (a_tg(1,isp1,ia1).gt.0.) then
            kg = log(2.) / a_tg(1,isp1,ia1)
        else
            kg = -9.
        endif

        a_nmXt = a_nmT(1,isp1,ia1) - a_nmS(1,isp1,ia1) - a_nmN(1,isp1,ia1)
        a_nmXSt = a_nmXt + a_nmS(1,isp1,ia1)
        
        ig1 = a_gid(1,isp1,ia1)
        a_GCTt = dble(a_ntnk(1,isp1,ia1,3)+a_ntnk(1,isp1,ia1,4))/g_ntn(ig1)
        a_GCXt = dble(a_ntnXk(1,isp1,ia1,3)+a_ntnXk(1,isp1,ia1,4))/g_ntnX(ig1)
        a_GCSt = dble(a_ntnSk(1,isp1,ia1,3)+a_ntnSk(1,isp1,ia1,4))/g_ntnS(ig1)
        a_GCXSt = dble(a_ntnXk(1,isp1,ia1,3)+a_ntnXk(1,isp1,ia1,4)+a_ntnSk(1,isp1,ia1,3)+a_ntnSk(1,isp1,ia1,4))/(g_ntnX(ig1)+g_ntnS(ig1))
        a_GCNt = dble(a_ntnNk(1,isp1,ia1,3)+a_ntnNk(1,isp1,ia1,4))/g_ntnN(ig1)

        write(271,271) t, &
            dble(ia1), &
            dble(a_sa(1,isp1,ia1)), &
            dble(a_sn(1,isp1,ia1)), &
            dble(a_rsn(1,isp1,ia1)), & 
            a_sr(1,isp1,ia1), &
            dble(a_nd(1,isp1,ia1)), &
            dble(a_nmt(1,isp1,ia1)), &
            dble(a_nmS(1,isp1,ia1)), &
            dble(a_nmN(1,isp1,ia1)), &
            dble(a_nrt(1,isp1,ia1)), &
            dble(a_nrS(1,isp1,ia1)), &
            dble(a_nrN(1,isp1,ia1)), &
            a_tb(1,isp1,ia1), &
            a_tg(1,isp1,ia1), &
            dble(isp1), &
            dble(a_gid(1,isp1,ia1)), &
            dble(a_cn(1,isp1,ia1)), &
            dble(a_ck(1,isp1,ia1,1)), &
            dble(a_cnt(1,isp1,ia1,1)), &
            dble(a_cXSN(1,isp1,ia1,1)), &
            (a_ntnk(1,isp1,ia1,nt1)/dble(dimnnt),nt1=1,4), &
            a_q0DNAC(1,isp1,ia1), &
            a_q0DNAN(1,isp1,ia1), &
            a_q0DNAP(1,isp1,ia1), &
            a_q0aaC(1,isp1,ia1), &
            a_q0aaN(1,isp1,ia1), &
            a_qC(1,isp1,ia1), &
            a_qN(1,isp1,ia1), &
            a_qP(1,isp1,ia1), &
            kg, &
            a_GCTt, &
            a_GCXt, &
            a_GCSt, &
            a_GCXSt, &
            a_GCNt, &
            dble(a_nmXt), &
            dble(a_nmXSt), &
            (dble(a_naa(1,isp1,ia1,iaa1)) / dble(g_naat(ig1)), iaa1 = 1, 20), &
            -9., &
            -9., &
            -9.
       
    endif
enddo ! ia1
enddo ! isp1

271 format(64e24.16)
    
close(271)

tpnPOP = tpnPOP + dtpPOP

endif

!    
! ------------------
! --- output ORF ---
! ------------------
!
if ((opORF.ge.0).and.(t.ge.tpnORF)) then

write(*,*) 'IAM\ORF: Doing protein stats.'

! note: this calculation does not account for fixed changes
do ig1 = 1, ngt
    do ip1 = 1, g_np(ig1)
        g_pnmS(ig1,ip1) = 0.
        g_pnmN(ig1,ip1) = 0.
        g_pnrS(ig1,ip1) = 0.
        g_pnrN(ig1,ip1) = 0.
    enddo
enddo

if ((oHANS.ge.1).and.(opORF.le.1)) then
    SRT = 0.
    naa = 0
    do isp1 = 1, dimnsp
    do ia1 = 1, nat_sp(1,isp1)
    if (a_sa(1,isp1,ia1).ge.1) then
        SRT = SRT + a_sr(1,isp1,ia1)
        naa = naa + 1
    endif
    enddo
    enddo
endif

if (opORF.eq.0) then
    npORF = naa
endif

SRTp = 0.
isp1 = 1
ia2 = 0
do ia1 = 1, npORF
7710 continue
    if (opORF.eq.1) then
        if (oHANS.ge.1) then
            r = r4_uni(jdum)
            SRCDF = 0.
            do isp1 = 1, dimnsp
            do ia2 = 1, nat_sp(1,isp1)
            if (a_sa(1,isp1,ia2).ge.1) then
                SRCDF = SRCDF + a_sr(1,isp1,ia2) / SRT
                if (SRCDF.ge.r) then
                    goto 7720
                endif
            endif
            enddo
            enddo
            write(*,*) 'IAM\ORF: ERR on agent selection.'
            goto 9100
7720 continue
        else
            r = r4_uni(jdum)
            isp1 = 1 + r * dimnsp
            r = r4_uni(jdum)
            ia2 = 1 + r * nat_sp(1,isp1)
        endif
    else
        ia2 = ia2 + 1
        if (ia2.gt.nat_sp(1,isp1)) then
            isp1 = isp1 + 1
            ia2 = 1
        endif
        if (isp1.gt.dimnsp) then
            write(*,*) 'IAM\ORF: Warning: Not enough active agents.'
            goto 7730
        endif
    endif
    ia3 = ia2
    !write(*,*) 'IAM\ORF: ia1, ia2, ia3 = ', ia1, ia2, ia3
    if (a_sa(1,isp1,ia3).eq.0) then
        !write(*,*) 'IAM\ORF: Inactive agent, trying again...'
        goto 7710
    endif

    SRTp = SRTp + a_sr(1,isp1,ia3)
    ig1 = a_gid(1,isp1,ia3)
    do ic1 = 1, a_cn(1,isp1,ia3)
        int1 = a_ck(ic1,isp1,ia3,1)
        do ip1 = 1, g_np(ig1)
            if ((int1.ge.g_pntstart(ig1,ip1)).and.(int1.le.g_pntstop(ig1,ip1))) then
                if (a_cmr(ic1,isp1,ia3,1).eq.1) then
                    if (a_cXSN(ic1,isp1,ia3,1).eq.1) then
                        g_pnmS(ig1,ip1) = g_pnmS(ig1,ip1) + a_sr(1,isp1,ia3)
                    elseif (a_cXSN(ic1,isp1,ia3,1).eq.2) then
                        g_pnmN(ig1,ip1) = g_pnmN(ig1,ip1) + a_sr(1,isp1,ia3)
                    else
                        write(*,*) 'IAM: ERR: Invalid change at protein (S/N) (1).'
                        write(*,*) 'IAM: isp1, ia3 = ', isp1, ia3
                        write(*,*) 'IAM: ic1, int1, cXSN = ', ic1, int1, a_cXSN(ic1,isp1,ia3,1)
                        write(*,*) 'IAM: ip1 = ', ip1
                        write(*,*) 'IAM: start, stop, other = ', g_pntstart(ig1,ip1), g_pntstop(ig1,ip1), g_pother(ig1,ip1)
                        goto 9100
                    endif
                elseif (a_cmr(ic1,isp1,ia3,1).eq.2) then
                    if (a_cXSN(ic1,isp1,ia3,1).eq.1) then
                        g_pnrS(ig1,ip1) = g_pnrS(ig1,ip1) + a_sr(1,isp1,ia3)
                    elseif (a_cXSN(ic1,isp1,ia3,1).eq.2) then
                        g_pnrN(ig1,ip1) = g_pnrN(ig1,ip1) + a_sr(1,isp1,ia3)
                    else
                        write(*,*) 'IAM: ERR: Invalid change at protein (S/N) (2).'
                        write(*,*) 'IAM: isp1, ia3 = ', isp1, ia3
                        write(*,*) 'IAM: ic1, int1, cXSN = ', ic1, int1, a_cXSN(ic1,isp1,ia3,1)
                        write(*,*) 'IAM: ip1 = ', ip1
                        write(*,*) 'IAM: start, stop, other = ', g_pntstart(ig1,ip1), g_pntstop(ig1,ip1), g_pother(ig1,ip1)
                        goto 9100
                    endif
                endif
            endif ! protein found
        enddo ! ip1
    enddo ! ic1

enddo ! ia1

7730 continue

!pnmta = 0.
pnmSa = 0.
pnmNa = 0.
!pnrta = 0.
pnrSa = 0.
pnrNa = 0.

ppiT = 0.
ppiSa = 0.
ppiNa = 0.
ppiMa = 0.
ppiRa = 0.
ppiTa = 0.

do ig1 = 1, ngt
    write(*,*) 'IAM: ig1 = ', ig1
    do ip1 = 1, g_np(ig1)
        g_pnmS(ig1,ip1) = g_pnmS(ig1,ip1) / SRTp
        g_pnmN(ig1,ip1) = g_pnmN(ig1,ip1) / SRTp
        g_pnrS(ig1,ip1) = g_pnrS(ig1,ip1) / SRTp
        g_pnrN(ig1,ip1) = g_pnrN(ig1,ip1)  / SRTp
        pnmSa = pnmSa + g_pnmS(ig1,ip1)
        pnmNa = pnmNa + g_pnmN(ig1,ip1)
        pnrSa = pnrSa + g_pnrS(ig1,ip1)
        pnrNa = pnrNa + g_pnrN(ig1,ip1)
        if (g_pntnS(ig1,ip1).gt.0.) then
            ppiS = (g_pnmS(ig1,ip1)+g_pnrS(ig1,ip1))/g_pntnS(ig1,ip1)
        else
            ppiS = -9.
        endif
        if (g_pntnN(ig1,ip1).gt.0.) then
            ppiN = (g_pnmN(ig1,ip1)+g_pnrN(ig1,ip1))/g_pntnN(ig1,ip1)
        else
            ppiN = -9.
        endif
        ppiM = (g_pnmS(ig1,ip1)+g_pnmN(ig1,ip1))/g_pntn(ig1,ip1)
        ppiR = (g_pnrS(ig1,ip1)+g_pnrN(ig1,ip1))/g_pntn(ig1,ip1)
        ppiT = ppiM+ppiR
        ppiSa = ppiSa + ppiS
        ppiNa = ppiNa + ppiN
        ppiMa = ppiMa + ppiM 
        ppiRa = ppiRa + ppiR 
        ppiTa = ppiTa + ppiT
    enddo

    pnmta = (pnmSa+pnmNa) / g_np(ig1)
    pnmSa = pnmSa / g_np(ig1)
    pnmNa = pnmNa / g_np(ig1)
    pnrta = (pnrSa+pnrNa) / g_np(ig1)
    pnrSa = pnrSa / g_np(ig1)
    pnrNa = pnrNa / g_np(ig1)
    ppiSa = ppiSa / g_np(ig1)
    ppiNa = ppiNa / g_np(ig1)
    ppiMa = ppiMa / g_np(ig1)
    ppiRa = ppiRa / g_np(ig1)
    ppiTa = ppiTa / g_np(ig1)

    write(*,*) 'IAM:'
    write(*,*) 'IAM: Protein statistics:'
    write(*,*) 'IAM: ppiS(low) = ', ppiSa
    write(*,*) 'IAM: ppiN(low) = ', ppiNa
    write(*,*) 'IAM: ppiM(low) = ', ppiMa
    write(*,*) 'IAM: ppiR(low) = ', ppiRa
    write(*,*) 'IAM: ppiT(low) = ', ppiTa

enddo

tpnORF = tpnORF + dtpORF

endif  
            
!
! ------------
! --- NOTE ---
! ------------
!
if (t.ge.tNoteNext) then
    write(*,*) 'IAM\NOTE1: t (d,y) = ', t, t/365.25
    write(*,*) 'IAM\NOTE2: nat_sp_mx, ncg_sp_mx = ', nat_sp_mx, ncg_sp_mx
    write(*,*) 'IAM\NOTE3: naa, SRT = ', naa, SRT
    write(*,*) 'IAM\NOTE4: cp_delta (h,h/ky) = ', cp_delta/3600., (cp_delta/3600.)/(t/365.25/1000.)
    if (fmutrecn.gt.0.) then
        write(*,*) 'IAM\NOTE5: fmut, frec = ', fmutt/fmutrecn, frect/fmutrecn
    else
        write(*,*) 'IAM\NOTE5: fmut, frec < fmutrecn = 0'
    endif
    write(*,*) 'IAM\NOTE6.1: nadt, nsrdt = ', nadt, nsrdt
    write(*,*) 'IAM\NOTE6.2: nawt, nsrwt = ', nawt, nsrwt
    write(*,*) 'IAM\NOTE6.3: namt, nsrmt = ', namt, nsrmt
    write(*,*) 'IAM\NOTE6.4: nart, nsrrt = ', nart, nsrrt
    call cp
    tNoteNext = t + dtNote
endif

cp_output_delta = cp_output_delta + omp_get_wtime() - cp_output_t1

!
! ----------------------
! --- sarkar presets ---
! ----------------------
!
if (isarkar.eq.0) then     
write(*,*) 'IAM\SARK: Populating preset arrays.'

isarkar = 1

kga = 0.
SRT = 0.

do isp1 = 1, dimnsp
!write(*,*) 'IAM\SARK: isp1 = ', isp1
do ia1 = 1, nat_sp(1,isp1)

if (a_sa(1,isp1,ia1).ge.1) then

    if (oLN.eq.1) then
        LNC = SC_sp(isp1) / (KmC + SC_sp(isp1))
        LNN = SN_sp(isp1) / (KmN + SN_sp(isp1))
        LNP = SP_sp(isp1) / (KmP + SP_sp(isp1))
    else
        LNC = 1.
        LNN = 1.
        LNP = 1.
    endif

    if (oSlim.eq.0) then
        call calc_V(isp1,ia1,LNC,LNN,LNP)
        call calc_q0(isp1,ia1,q0C,q0N,q0P)
        call calc_A(isp1,ia1,q0C,q0N,q0P,kg)
    else
        call calc_A_Slim(isp1,ia1,LNC,LNN,LNP,kg)
    endif

    !write(*,*) 'IAM\SARK: ia1, kg = ', ia1, kg

    if (kg.gt.0.) then
        kga = kga + kg * a_sr(1,isp1,ia1)
        SRT = SRT + a_sr(1,isp1,ia1)
    endif

endif
enddo
enddo

kg = kga / SRT
    
isp1 = 1
ia1 = 1

write(*,*) 'IAM\SARK: kg = ', kg

!sarkt = 0.
!sarkn = 0.
sarksr0 = 100. ! make input
xsr = exp(kg*dt)
do ia2 = 1, dimns
    if (a_uEg(1,isp1,ia1)*fM.gt.0.) then
        sarkarm(ia2) = sarkar(sarksr0*xsr*.5*a_uEg(1,isp1,ia1)*fM,sarksr0*xsr,idum_sp(isp1))
    else
        sarkarm(ia2) = 0.
    endif
    !sarkt = sarkt + sarkarm(ia2)
    !sarkn = sarkn + 1.
    if (uREg*fR.gt.0.) then
        sarkarr(ia2) = sarkar(sarksr0*xsr*.5*uREg*fR,sarksr0*xsr,idum_sp(isp1))
    else
        sarkarr(ia2) = 0.
    endif
enddo

!write(*,*) 'IAM\SARK: sarkave/sarksr0 = ', sarkt/sarkn/sarksr0

write(*,*) 'IAM\SARK: Done.'

endif

!
! --- subpopulation loop ---
!
cp_partot_t1 = omp_get_wtime()
!$OMP parallel num_threads(ntd) default(shared) PRIVATE(isp1,isp2, itx, ia1, ia2, ia3, ia4, int1, int2, int3, int4, int5, ig1, ig2, r, PD, &
 P1, SF, LNC, LNN, LNP, mqC, mqN, mqP, q0C, q0N, q0P, kmut, krec, qmaxC, qmaxN, qmaxP, dtad, SRCDFd, SRCDFs, SRCDFr, SRCDF, &
 nt1, nt2, ic1, ic2, ic3, icn1, icn2, iXSN, irtype, ip1, nta, ntb, nta_o, nt5, pntstartx, pntstopx, pntstepx, cod1, cod2, cod1_o, &
iaa1, iaa2, iaa1_o, Loth, fpi, si, rni, lR, Rstart, Rend, kgave, PS, SUMVmax,kgC, kgN, kgP, kg, sr0, xsr, qaveC, qaveN, qaveP, &
nx, hx, rx, ix, fx, ierr, ncg_sp_prev,ixmr,srx,srmx,srrx,sax,xMuse,xRuse)
!$OMP DO
do isp1 = 1, dimnsp
!write(*,*) 'IAM: Subpopulation start: isp1 = ', isp1
if (isp1.ne.omp_get_thread_num()+1) then
    write(*,*) 'IAM: OMP error.'
    write(*,*) 'IAM: isp1, thread = ', isp1, omp_get_thread_num()
endif

!
! -------------------------
! --- adopt master time ---
! -------------------------
!
t_sp(isp1) = t
tDNext_sp(isp1) = tDNext

!
! ------------------------
! --- adopt master gsn ---
! ------------------------
!
gsn_sp(isp1) = gsn

!
! ----------------------
! --- reset counters ---
! ----------------------
!
ncg_sp(isp1) = 0
nadt_sp(isp1) = 0.
nsrdt_sp(isp1) = 0.
fmutt_sp(isp1) = 0.
frect_sp(isp1) = 0.
fmutrecn_sp(isp1) = 0.
namt_sp(isp1) = 0.
nsrmt_sp(isp1) = 0.
nmTt_sp(isp1) = 0.
nmSt_sp(isp1) = 0.
nmNt_sp(isp1) = 0.
nmXt_sp(isp1) = 0.
do nt1 = 1, 4
    do nt2 = 1, 4
        nmTij_sp(isp1,nt1,nt2) = 0.
        nmSij_sp(isp1,nt1,nt2) = 0.
        nmNij_sp(isp1,nt1,nt2) = 0.
        nmXij_sp(isp1,nt1,nt2) = 0.
    enddo
enddo
nart_sp(isp1) = 0.
nsrrt_sp(isp1) = 0.
nrTt_sp(isp1) = 0.
nrSt_sp(isp1) = 0.
nrNt_sp(isp1) = 0.
nrXt_sp(isp1) = 0.
do nt1 = 1, 4
    do nt2 = 1, 4
        nrTij_sp(isp1,nt1,nt2) = 0.
        nrSij_sp(isp1,nt1,nt2) = 0.
        nrNij_sp(isp1,nt1,nt2) = 0.
        nrXij_sp(isp1,nt1,nt2) = 0.
    enddo
enddo
rdeltat_sp(isp1) = 0.
do int1 = 1, 9
    rtype_sp(isp1,int1) = 0.
enddo
rnut_sp(isp1) = 0.
rn_sp(isp1) = 0.
fpt_sp(isp1) = 0.
fpn2_sp(isp1) = 0.
nawt_sp(isp1) = 0.
nsrwt_sp(isp1) = 0.

!
! -----------------------
! --- inner time loop ---
! -----------------------
!
do itx = 1, ntx
!write(*,*) 'IAM: Inner time start: isp1, itx, ntx = ', isp1, itx, ntx

!
! ------------------------------
! --- init extracellular der ---
! ------------------------------
!
RSC_sp(isp1) = 0.
RSN_sp(isp1) = 0.
RSP_sp(isp1) = 0.

!
! --------------------------
! --- inflow and outflow ---
! --------------------------
!
RSC_sp(isp1) = RSC_sp(isp1) + (Q_sp(isp1) * SinC - Q_sp(isp1) * SC_sp(isp1)) / VOL_sp(isp1)
RSN_sp(isp1) = RSN_sp(isp1) + (Q_sp(isp1) * SinN - Q_sp(isp1) * SN_sp(isp1)) / VOL_sp(isp1)
RSP_sp(isp1) = RSP_sp(isp1) + (Q_sp(isp1) * SinP - Q_sp(isp1) * SP_sp(isp1)) / VOL_sp(isp1)

!
! ---------------------------
! --- nutrient limitation ---
! ---------------------------
!
if (oLN.eq.1) then
    LNC = SC_sp(isp1) / (KmC + SC_sp(isp1))
    LNN = SN_sp(isp1) / (KmN + SN_sp(isp1))
    LNP = SP_sp(isp1) / (KmP + SP_sp(isp1))
else
    LNC = 1.
    LNN = 1.
    LNP = 1.
endif

!
! -----------------------
! --- main agent loop ---
! -----------------------
!
nar_sp(isp1) = 0

VmaxCa_sp(isp1) = 0.
VmaxNa_sp(isp1) = 0.
VmaxPa_sp(isp1) = 0.

if (oHANS.eq.-1) then ! oHANS

do ia1 = 1, nat_sp(1,isp1)
if (a_sa(1,isp1,ia1).ge.1) then
    a_srm(1,isp1,ia1) = 1.
    a_srr(1,isp1,ia1) = 1.
    nar_sp(isp1) = nar_sp(isp1) + 1
endif ! sa
enddo ! ia1
    
elseif (oHANS.eq.0) then ! oHANS

cp_ENS_t1_sp(isp1) = omp_get_wtime()
do ia1 = 1, nat_sp(1,isp1)
if (a_sa(1,isp1,ia1).ge.1) then

    call calc_V(isp1,ia1,LNC,LNN,LNP)
    call calc_q0(isp1,ia1,q0C,q0N,q0P)

    qmaxC = fqmaxC * q0C
    qmaxN = fqmaxN * q0N
    qmaxP = fqmaxP * q0P

    if (a_qC(1,isp1,ia1).ge.qmaxC) then
        a_VC(1,isp1,ia1) = 0.
    endif
    if (a_qN(1,isp1,ia1).ge.qmaxN) then
        a_VN(1,isp1,ia1) = 0.
    endif
    if (a_qP(1,isp1,ia1).ge.qmaxP) then
        a_VP(1,isp1,ia1) = 0.
    endif
    
    RSC_sp(isp1) = RSC_sp(isp1) - a_VC(1,isp1,ia1) / VOL_sp(isp1) * a_sr(1,isp1,ia1)
    RSN_sp(isp1) = RSN_sp(isp1) - a_VN(1,isp1,ia1) / VOL_sp(isp1) * a_sr(1,isp1,ia1)
    RSP_sp(isp1) = RSP_sp(isp1) - a_VP(1,isp1,ia1) / VOL_sp(isp1) * a_sr(1,isp1,ia1)

    ! ave Vmax for implicit    
    if (oImpExt.eq.1) then
        VmaxCa_sp(isp1) = VmaxCa_sp(isp1) + a_VmaxC(1,isp1,ia1) * a_sr(1,isp1,ia1)
        VmaxNa_sp(isp1) = VmaxNa_sp(isp1) + a_VmaxN(1,isp1,ia1) * a_sr(1,isp1,ia1)
        VmaxPa_sp(isp1) = VmaxPa_sp(isp1) + a_VmaxP(1,isp1,ia1) * a_sr(1,isp1,ia1)
    endif

    !
    ! quotas
    !
    a_qC(1,isp1,ia1) = a_qC(1,isp1,ia1) + a_VC(1,isp1,ia1) * dt
    a_qN(1,isp1,ia1) = a_qN(1,isp1,ia1) + a_VN(1,isp1,ia1) * dt
    a_qP(1,isp1,ia1) = a_qP(1,isp1,ia1) + a_VP(1,isp1,ia1) * dt

    !
    ! --- death & outflow ---
    !
    PD = (kd + Q_sp(isp1)/VOL_sp(isp1)) * dt
    if (PD.gt.0.) then
        r = r4_uni(idum_sp(isp1))
        if (r.le.PD) then
            !write(*,*) 'IAM: Death/outflow. ia1, sn = ', ia1, a_sn(1,isp1,ia1)
            nawt_sp(isp1) = nawt_sp(isp1) + 1.
            nsrwt_sp(isp1) = nsrwt_sp(isp1) + 1.
            call kill_agent(isp1,ia1)
            goto 399
        endif
    endif

    !
    ! --- division ---
    !
    if ((a_qC(1,isp1,ia1).ge.(2.*q0C)).and.(a_qN(1,isp1,ia1).ge.(2.*q0N)).and.(a_qP(1,isp1,ia1).ge.(2.*q0P))) then
        !write(*,*) 'IAM: Division. ia1, sn = ', ia1, a_sn(1,isp1,ia1)
        a_CLimit(1,isp1,ia1) = a_qC(1,isp1,ia1)/(2.*q0C)
        a_NLimit(1,isp1,ia1) = a_qN(1,isp1,ia1)/(2.*q0N)
        a_PLimit(1,isp1,ia1) = a_qP(1,isp1,ia1)/(2.*q0P)
        nadt_sp(isp1) = nadt_sp(isp1) + 1.
        nsrdt_sp(isp1) = nsrdt_sp(isp1) + 1.
        a_nd(1,isp1,ia1) = a_nd(1,isp1,ia1) + 1
        a_tg(1,isp1,ia1) = t_sp(isp1) - a_tb(1,isp1,ia1)
        a_tb(1,isp1,ia1) = t_sp(isp1)
        !correction for uptake at division
        if (oCAD.eq.1) then
            dtad = 0.
            if ((a_CLimit(1,isp1,ia1).le.a_NLimit(1,isp1,ia1)).and.(a_CLimit(1,isp1,ia1).le.a_PLimit(1,isp1,ia1)) &
                .and.(a_VC(1,isp1,ia1).gt.0.)) then
                dtad = (a_qC(1,isp1,ia1) - 2. * q0C) / a_VC(1,isp1,ia1)
            elseif ((a_NLimit(1,isp1,ia1).le.a_CLimit(1,isp1,ia1)).and.(a_NLimit(1,isp1,ia1).le.a_PLimit(1,isp1,ia1)) &
                .and.(a_VN(1,isp1,ia1).gt.0.)) then
                dtad = (a_qN(1,isp1,ia1) - 2. * q0N) / a_VN(1,isp1,ia1)
            elseif ((a_PLimit(1,isp1,ia1).le.a_CLimit(1,isp1,ia1)).and.(a_PLimit(1,isp1,ia1).le.a_NLimit(1,isp1,ia1)) &
                .and.(a_VP(1,isp1,ia1).gt.0.)) then
                dtad = (a_qP(1,isp1,ia1) - 2. * q0P) / a_VP(1,isp1,ia1)
            endif
            RSC_sp(isp1) = RSC_sp(isp1) - a_VC(1,isp1,ia1) / VOL_sp(isp1) * a_sr(1,isp1,ia1) * dtad / dt
            RSN_sp(isp1) = RSN_sp(isp1) - a_VN(1,isp1,ia1) / VOL_sp(isp1) * a_sr(1,isp1,ia1) * dtad / dt
            RSP_sp(isp1) = RSP_sp(isp1) - a_VP(1,isp1,ia1) / VOL_sp(isp1) * a_sr(1,isp1,ia1) * dtad / dt
            a_qC(1,isp1,ia1) = a_qC(1,isp1,ia1) + a_VC(1,isp1,ia1) * dtad
            a_qN(1,isp1,ia1) = a_qN(1,isp1,ia1) + a_VN(1,isp1,ia1) * dtad
            a_qP(1,isp1,ia1) = a_qP(1,isp1,ia1) + a_VP(1,isp1,ia1) * dtad
            a_tg(1,isp1,ia1) = a_tg(1,isp1,ia1) - dtad
            a_tb(1,isp1,ia1) = a_tb(1,isp1,ia1) - dtad
        endif

        SF = 0.5 * kreft_p(RCVmm, RCVmcv, RCVmcvl, RCVmll, RCVmhl, idum_sp(isp1),kn,fn,wn)
        mqC = a_qC(1,isp1,ia1)
        mqN = a_qN(1,isp1,ia1)
        mqP = a_qP(1,isp1,ia1)
        a_qC(1,isp1,ia1) = mqC * SF
        a_qN(1,isp1,ia1) = mqN * SF
        a_qP(1,isp1,ia1) = mqP * SF
       
        kmut = 0
        if (pcheck(a_uEg(1,isp1,ia1)*fM, idum_sp(isp1)).eq.1)then ! mutation
            !write(*,*) 'IAM\MUT: t_sp, ia1 = ', t_sp(isp1), ia1
            namt_sp(isp1) = namt_sp(isp1) + 1.
            nsrmt_sp(isp1) = nsrmt_sp(isp1) + 1.
            r = r4_uni(idum_sp(isp1)) ! decide if mother or daughter
            if (r.lt.0.5) then
                kmut = 1
            else
                kmut = 2
            endif
            if (kmut.eq.1) then
                a_srm(1,isp1,ia1) = 1.
            endif
        endif

        krec = 0
        if (pcheck(uREg*fR, idum_sp(isp1)).eq.1)then ! recombination
            !write(*,*) 'IAM\REC1: t_sp, ia1 = ', t_sp(isp1), ia1
            nart_sp(isp1) = nart_sp(isp1) + 1.
            nsrrt_sp(isp1) = nsrrt_sp(isp1) + 1.
            r = r4_uni(idum_sp(isp1)) ! decide if mother or daughter
            if (r.lt.0.5) then
                krec = 1
            else
                krec = 2
            endif
            if (krec.eq.1) then
                a_srr(1,isp1,ia1) = 1.
            endif
        endif
        
        if (oKillDaughter.eq.0) then

            call get_new_ia(isp1,ia2,ierr)
            if (ierr.eq.1) then
                write(*,*) 'IAM: ERR: nat_sp > dimna_sp.', nat_sp(1,isp1), dimna_sp
                err_sp(isp1) = 301
                goto 9050
            endif
            
            call copy_agent(1,isp1,ia1,1,isp1,ia2)
            
            a_qC(1,isp1,ia2) = mqC * (1. - SF)
            a_qN(1,isp1,ia2) = mqN * (1. - SF)
            a_qP(1,isp1,ia2) = mqP * (1. - SF)
            if (kmut.eq.2) then
                a_srm(1,isp1,ia2) = 1.
            else
                a_srm(1,isp1,ia2) = 0.
            endif
            if (krec.eq.2) then
                a_srr(1,isp1,ia2) = 1.
            else
                a_srr(1,isp1,ia2) = 0.
            endif
            
            if (a_srr(1,isp1,ia2).gt.0.) then ! should be or ia1?
                nar_sp(isp1) = nar_sp(isp1) + 1
            endif
        
            if (oSlim.eq.0) then
                call new_sn(1,isp1,ia2)
            endif

        endif

    endif ! division

399 continue
   
endif ! sa
enddo ! ia1
cp_ENS_delta_sp(isp1) = cp_ENS_delta_sp(isp1) + omp_get_wtime() - cp_ENS_t1_sp(isp1)

elseif (oHANS.eq.1) then ! oHANS
    
if (opD.eq.2) then
    goto 6310
endif

!
! --- growth ---
!
6210 continue
    
cp_HANS1_t1_sp(isp1) = omp_get_wtime()
SRT_sp(isp1) = 0.
do ia1 = 1, nat_sp(1,isp1)
if (a_sa(1,isp1,ia1).ge.1) then
    
    a_sa(1,isp1,ia1) = 1

    if (oSlim.eq.0) then
        call calc_V(isp1,ia1,LNC,LNN,LNP)
        call calc_q0(isp1,ia1,q0C,q0N,q0P)
        call calc_A(isp1,ia1,q0C,q0N,q0P,kg)
        if (kg.gt.0.) then
            a_tg(1,isp1,ia1) = log(2.) / kg
        else
            a_tg(1,isp1,ia1) = 1.e9
        endif
        a_tb(1,isp1,ia1) = t_sp(isp1)
    else
        call calc_A_Slim(isp1,ia1,LNC,LNN,LNP,kg)
    endif
    
    sr0 = a_sr(1,isp1,ia1)
    SRT_sp(isp1) = SRT_sp(isp1) + sr0
    a_srm(1,isp1,ia1) = 0.
    a_srr(1,isp1,ia1) = 0.

    if (kg.gt.0.) then

        xsr = exp(kg*dt)
        a_sr(1,isp1,ia1) = sr0 * xsr
        SRT_sp(isp1) = SRT_sp(isp1) + (a_sr(1,isp1,ia1) - sr0)
        nsrdt_sp(isp1) = nsrdt_sp(isp1) + (a_sr(1,isp1,ia1) - sr0)
        a_nd(1,isp1,ia1) = a_nd(1,isp1,ia1) + (a_sr(1,isp1,ia1) - sr0)
        
        if (oHANSs.eq.1) then
            nx = a_sr(1,isp1,ia1)
            hx = dt * kg / log(2.) * .5 * a_uEg(1,isp1,ia1)*fM * a_sr(1,isp1,ia1)
            rx = bucci(hx, nx, idum_sp(isp1))
            a_srm(1,isp1,ia1) = rx
            hx = dt * kg / log(2.) * .5 * uREg*fR * a_sr(1,isp1,ia1)
            rx = bucci(hx, nx, idum_sp(isp1))
            a_srr(1,isp1,ia1) = rx
        elseif (oHANSs.eq.2) then
            if (a_uEg(1,isp1,ia1)*fM.gt.0.) then
                a_srm(1,isp1,ia1) = sarkar(sr0*xsr*.5*a_uEg(1,isp1,ia1)*fM,sr0*xsr,idum_sp(isp1))
            else
                a_srm(1,isp1,ia1) = 0.
            endif
            if (uREg*fR.gt.0.) then
                a_srr(1,isp1,ia1) = sarkar(sr0*xsr*.5*uREg*fR,sr0*xsr,idum_sp(isp1))
            else
                a_srr(1,isp1,ia1) = 0.
            endif
        else
            if (a_uEg(1,isp1,ia1)*fM.gt.0.) then
                r = r4_uni(idum_sp(isp1))
                a_srm(1,isp1,ia1) = sr0/sarksr0*dble(sarkarm(min(1+int(r*dimns),dimns)))
            else
                a_srm(1,isp1,ia1) = 0.
            endif
            if (uREg*fR.gt.0.) then
                r = r4_uni(idum_sp(isp1))
                a_srr(1,isp1,ia1) = sr0/sarksr0*dble(sarkarr(min(1+int(r*dimns),dimns)))
            else
                a_srr(1,isp1,ia1) = 0.
            endif
        endif

        !mutation/recombination
        !events:
        !nsrmt_sp(isp1) = nsrmt_sp(isp1) + (a_sr(1,isp1,ia1) - sr0) * a_uEg(1,isp1,ia1)*fM
        !nsrrt_sp(isp1) = nsrrt_sp(isp1) + (a_sr(1,isp1,ia1) - sr0) * uREg*fR
        !cells:
        nsrmt_sp(isp1) = nsrmt_sp(isp1) + a_srm(1,isp1,ia1)
        nsrrt_sp(isp1) = nsrrt_sp(isp1) + a_srr(1,isp1,ia1)

    endif ! kg>0

    fmutt_sp(isp1) = fmutt_sp(isp1) + a_srm(1,isp1,ia1)/a_sr(1,isp1,ia1)
    frect_sp(isp1) = frect_sp(isp1) + a_srr(1,isp1,ia1)/a_sr(1,isp1,ia1)
    fmutrecn_sp(isp1) = fmutrecn_sp(isp1) + 1.
    
endif ! sa
enddo ! ia1

cp_HANS1_delta_sp(isp1) = cp_HANS1_delta_sp(isp1) + omp_get_wtime() - cp_HANS1_t1_sp(isp1)

if (opD.eq.2) then
    goto 6410
endif

!
! --- dilution & mut and reco agents ---
!
6310 continue
cp_HANS2_t1_sp(isp1) = omp_get_wtime()

PS = Popic_sp(isp1) / SRT_sp(isp1)

do ia1 = 1, nat_sp(1,isp1)
if (a_sa(1,isp1,ia1).eq.1) then

    !write(*,*) 'IAM\DIL: Start isp1, ia1 = ', isp1, ia1
    
    srx = a_sr(1,isp1,ia1)
    srmx = a_srm(1,isp1,ia1)
    srrx = a_srr(1,isp1,ia1)
    sax = 1
    a_sr(1,isp1,ia1) = 0.
    a_srm(1,isp1,ia1) = 0.
    a_srr(1,isp1,ia1) = 0.

    !cp_HANS2a_t1_sp(isp1) = omp_get_wtime()

    !
    ! wildtype
    !
    nx = srx
    hx = PS * nx
    if (PS.ge.1.) then
        rx = nx
    else
        rx = bucci(hx, nx, idum_sp(isp1))
    endif
    fx = rx / nx !for proportional, see below
    ix = nint(rx)
    if (ix.eq.0) then
        nawt_sp(isp1) = nawt_sp(isp1) + 1.
    endif
    if (oHANSr.eq.1) then
        srx = dble(ix)
    else
        srx = rx
    endif
    nsrwt_sp(isp1) = nsrwt_sp(isp1) + (nx-rx)

    !cp_HANS2a_delta_sp(isp1) = cp_HANS2a_delta_sp(isp1) + omp_get_wtime() - cp_HANS2a_t1_sp(isp1)
    
    !
    ! mutations
    !
    !cp_HANS2b_t1_sp(isp1) = omp_get_wtime()    
    if (oHANSd.eq.0) then
        nx = srmx
        hx = PS * nx
        if (PS.ge.1.) then
            rx = nx
        else
            rx = bucci(hx, nx, idum_sp(isp1))
        endif
    else
        rx = srmx * fx
    endif

    ix = nint(rx)
    if (ix.ge.1) then
        srx = srx - dble(ix)
        do ia3 = 1, ix
            
            if ((sax.eq.1).and.(srx.lt.1.)) then ! recycle mother if empty
                ia2 = ia1
                sax = 0
                !write(*,*) 'IAM\DIL: Recycle mother.'
            else
                call get_new_ia(isp1,ia2,ierr)
                if (ierr.eq.1) then
                    write(*,*) 'IAM: ERR: nat_sp > dimna_sp.', nat_sp(1,isp1), dimna_sp
                    err_sp(isp1) = 401
                    goto 9050
                endif
                call copy_agent(1,isp1,ia1,1,isp1,ia2)
                !write(*,*) 'IAM\DIL: New agent.'
            endif

            a_sa(1,isp1,ia2) = 2
            a_srm(1,isp1,ia2) = 1.
            a_srr(1,isp1,ia2) = 0.

            if (oSlim.eq.0) then
                call new_sn(1,isp1,ia2)
            endif

            namt_sp(isp1) = namt_sp(isp1) + 1.
            
            if (oHANSu.eq.1) then
                a_sr(1,isp1,ia2) = dble(ix)
                goto 6315
            else
                a_sr(1,isp1,ia2) = 1.
            endif
            
        enddo ! ia3

    endif ! ix >= 1

6315 continue

    !cp_HANS2b_delta_sp(isp1) = cp_HANS2b_delta_sp(isp1) + omp_get_wtime() - cp_HANS2b_t1_sp(isp1)
    
    !
    ! recombinations
    !
    !cp_HANS2c_t1_sp(isp1) = omp_get_wtime()
    
    if (oHANSd.eq.0) then
        nx = srrx
        hx = PS * nx
        if (PS.ge.1.) then
            rx = nx
        else
            rx = bucci(hx, nx, idum_sp(isp1))
        endif
    else
        rx = srrx * fx
    endif

    ix = nint(rx)
    if (ix.ge.1) then

        srx = srx - dble(ix)

        do ia3 = 1, ix
            
            if ((sax.eq.1).and.(srx.lt.1.)) then ! recycle mother if empty
                ia2 = ia1
                sax = 0
                !write(*,*) 'IAM\DIL: Recycle mother.'
            else
                call get_new_ia(isp1,ia2,ierr)
                if (ierr.eq.1) then
                    write(*,*) 'IAM: ERR: nat_sp > dimna_sp.', nat_sp(1,isp1), dimna_sp
                    err_sp(isp1) = 501
                    goto 9050
                endif
                call copy_agent(1,isp1,ia1,1,isp1,ia2)
                !write(*,*) 'IAM\DIL: New agent.'
            endif

            a_sa(1,isp1,ia2) = 2
            a_srm(1,isp1,ia2) = 0.
            a_srr(1,isp1,ia2) = 1.

            if (oSlim.eq.0) then
                call new_sn(1,isp1,ia2)
            endif

            nar_sp(isp1) = nar_sp(isp1) + 1
            nart_sp(isp1) = nart_sp(isp1) + 1.
            
            if (oHANSu.eq.1) then
                a_sr(1,isp1,ia2) = dble(ix)
                goto 6325
            else
                a_sr(1,isp1,ia2) = 1.
            endif
            
        enddo

    endif

6325 continue

    if (sax.eq.1) then  ! if mother not recycled
        if (srx.lt.1.) then  ! if mother empty
            call kill_agent(isp1,ia1)
        else
            a_sr(1,isp1,ia1) = srx
        endif
    endif
     
    !write(*,*) 'IAM\DIL: Stop isp1, ia1 = ', isp1, ia1
    
endif ! sa
enddo ! ia1

!write(*,*) 'IAM\DIL: Done.'

cp_HANS2_delta_sp(isp1) = cp_HANS2_delta_sp(isp1) + omp_get_wtime() - cp_HANS2_t1_sp(isp1)
if (opD.eq.2) then
    goto 6210
endif

6410 continue
     
endif

!
! -------------------------
! --- update cell count ---
! -------------------------
!
SRT_sp(isp1) = 0.
if ((nat_sp(1,isp1)-naf_sp(isp1)).eq.0) then
    write(*,*) 'IAM\UCC: Warning: Empty sp. isp1 = ', isp1
    write(*,*) 'IAM\UCC: isp1, nat_sp = ', isp1, nat_sp(1,isp1)
    write(*,*) 'IAM\UCC: isp1, SRT_sp = ', isp1, SRT_sp(isp1)
    goto 9050
endif
do ia1 = 1, nat_sp(1,isp1)
    if (a_sa(1,isp1,ia1).ge.1) then
        SRT_sp(isp1) = SRT_sp(isp1) + a_sr(1,isp1,ia1)
    endif
enddo

!
! ----------------------------------
! --- mutations & recombinations ---
! ----------------------------------
!
! --- donor selection for recombination ---
!
cp_mr1_t1_sp(isp1) = omp_get_wtime()
if (nar_sp(isp1).gt.0) then
    if (oHANS.ge.1) then
        r = r4_uni(idum_sp(isp1))
        SRCDFd = 1. / dble(nar_sp(isp1)*xR)
        SRCDFs = r * SRCDFd
        SRCDFr = SRCDFs
        SRCDF = 0.
        ia1 = 1
        do ia2 = 1, nat_sp(1,isp1)
            if (a_sa(1,isp1,ia2).ge.1) then
                SRCDF = SRCDF + a_sr(1,isp1,ia2) / SRT_sp(isp1)
701 continue
                if (SRCDF.ge.SRCDFr) then
                    iad_sp(isp1,ia1) = ia2
                    if (ia1.eq.(nar_sp(isp1)*xR)) then
                        goto 702
                    endif
                    ia1 = ia1 + 1
                    SRCDFr = SRCDFr + SRCDFd
                    goto 701
                endif
            endif
        enddo
        iad_sp(isp1,ia1) = iad_sp(isp1,ia1-1) ! last one may not catch due to rounding error
        if (ia1.ne.(nar_sp(isp1)*xR)) then
            write(*,*) 'IAM\REC: ERR selecting donor.'
            write(*,*) 'IAM\REC: isp1 = ', isp1
            write(*,*) 'IAM\REC: ia1, nar_sp = ', ia1, nar_sp(isp1)
            err_sp(isp1) = 801
            goto 9050
        endif
702 continue
    endif ! oHANS
endif ! nar>0
cp_mr1_delta_sp(isp1) = cp_mr1_delta_sp(isp1) + omp_get_wtime() - cp_mr1_t1_sp(isp1)

!
! --- agent loop ---
!
cp_mr2_t1_sp(isp1) = omp_get_wtime()
do ia1 = 1, nat_sp(1,isp1)
if (a_sa(1,isp1,ia1).ge.1) then

ncg_sp_prev = ncg_sp(isp1)

!
! --- mutation ---
!
if (a_srm(1,isp1,ia1).ge.1.) then

    !if (isp1.eq.1) then
    !if (a_sn(1,isp1,ia1).eq.390) then
    !write(*,*) 'IAM\MUT: t_sp, ia1 = ', t_sp(isp1), ia1
    !endif

    ! kill mutant
    if (fmKill.gt.0.) then
        r = r4_uni(idum_sp(isp1))
        if (r.lt.fmKill) then
            call kill_agent(isp1,ia1)
            goto 403
        endif
    endif
    
    ig1 = a_gid(1,isp1,ia1)
    a_srm(1,isp1,ia1) = 0.

    !if (oHANS.eq.-1) then
    !    xMuse = bucci(dble(xM), dble(xM*10), idum_sp(isp1))
    !else
        xMuse = xM
    !endif
    
    do ixmr = 1, xMuse
    !if (isp1.eq.1) then
    !if (a_sn(1,isp1,ia1).eq.390) then
    !write(*,*) 'IAM\MUT: ixmr, xMuse = ', ixmr, xMuse
    !endif
    
    !
    ! select from class
    !
    r = r4_uni(idum_sp(isp1))
    P1 = 0. 
    do int2 = 1, 3
        P1 = P1 + a_ntnk(1,isp1,ia1,int2) * uEbk(int2) / a_uEg(1,isp1,ia1)
        if (P1.ge.r) then
            goto 302
        endif
    enddo
302 continue
    !if (isp1.eq.1) then
    !if (a_sn(1,isp1,ia1).eq.390) then
    !write(*,*) 'IAM\MUT: P1 = ', P1
    !endif

    !
    ! select from bp
    !
    int3 = 0
303 continue
    r = r4_uni(idum_sp(isp1))
    int1 = min(int(1+r*dimnnt),dimnnt)
    if (g_nt(ig1,int1).ne.int2) then
        int3 = int3 + 1
        if (int3.gt.dimnnt) then ! if dimnnt random tries are unsucessful
            write(*,*) 'IAM\MUT: Warning: Going to sequential bp selection.'
            do int4 = 1, dimnnt ! go sequentially starting from last try position
                int1 = int1 + 1
                if (int1.gt.dimnnt) then
                    int1 = 1
                endif
                if (g_nt(ig1,int1).eq.int2) then
                    goto 304
                endif
            enddo
            write(*,*) 'IAM\MUT: ERR selecting bp.'
            err_sp(isp1) = 901
            goto 9050
        endif
        goto 303
    endif
304 continue
    !if (isp1.eq.1) then
    !if (a_sn(1,isp1,ia1).eq.390) then
    !write(*,*) 'IAM\MUT: Changing int1 ', int1
    !endif
    ! get context
    do int5 = int1 - 2, int1 + 2
        int3 = int5 - int1 + 3
        nt5(int3) = -9
        if ((int5.ge.1).and.(int5.le.dimnnt)) then
            nt5(int3) = g_nt(ig1,int5)
        endif
    enddo
    ! update for changes
    icn1 = -1
    icn2 = 1
    do ic1 = 1, a_cn(1,isp1,ia1)
        int5 = a_ck(ic1,isp1,ia1,1)
        if ((int5.ge.(int1-2)).and.(int5.le.(int1+2))) then
            int3 = int5 - int1 + 3
            nt5(int3) = a_cnt(ic1,isp1,ia1,1)
            if (int5.eq.int1) then
                icn1 = ic1
                icn2 = 0
                !if (isp1.eq.1) then
                !if (a_sn(1,isp1,ia1).eq.390) then
                !write(*,*) 'IAM\MUT: Existing change at this location (1).'
                !write(*,*) 'IAM\MUT: isp1, ia1 = ', isp1, ia1
                !write(*,*) 'IAM\MUT: ic1, int1 = ', ic1, int1
                !endif
                if (oHANS.eq.-1) then
                    goto 303
                endif
            endif
        endif
    enddo
    nt1 = nt5(3)
    !if (isp1.eq.1) then
    !if (a_sn(1,isp1,ia1).eq.390) then
    !write(*,*) 'IAM\MUT: From ', nt1, int_snt(nt1)
    !endif
    !
    ! select destination
    !
    r = r4_uni(idum_sp(isp1))
    P1 = 0.
    do int2 = 1, 2
        nt2 = cint(nt1,int2)
        P1 = P1 + uEbij(nt1,nt2) / uEbk(nt1)
        if (P1.ge.r) then
            goto 306
        endif
    enddo
    nt2 = cint(nt1,int2)
306 continue
    !if (isp1.eq.1) then
    !if (a_sn(1,isp1,ia1).eq.390) then
    !write(*,*) 'IAM\MUT: To ', nt2, int_snt(nt2)
    !endif
    !
    !
    ! record change
    !
    call rec_c(isp1,ia1,int1,nt5,nt2,1,icn1,iXSN,ierr,-9.,-9.)
    if (ierr.eq.1) then
        write(*,*) 'IAM\MUT: ERR 1 in rec_c: jcn > dimnc.', a_cn(1,isp1,ia1), dimnc
        write(*,*) 'IAM\MUT: isp1, ia1 = ', isp1, ia1
        !err_sp(isp1) = 1001
        !goto 9050
        call kill_agent(isp1,ia1)
        ncg_sp(isp1) = ncg_sp_prev
        goto 403
    endif
    if (ierr.eq.2) then
        write(*,*) 'IAM\MUT: ERR 2 in rec_c. Invalid iNS.'
        err_sp(isp1) = 1101
        goto 9050
    endif
    if (a_sa(1,isp1,ia1).lt.1) then
        namt_sp(isp1) = namt_sp(isp1) - 1
        nsrmt_sp(isp1) = nsrmt_sp(isp1) - a_sr(1,isp1,ia1)
        nawt_sp(isp1) = nawt_sp(isp1) + 1
        nsrwt_sp(isp1) = nsrwt_sp(isp1) + a_sr(1,isp1,ia1)
        goto 403
    endif
    !if (isp1.eq.1) then
    !if (a_sn(1,isp1,ia1).eq.390) then
    !write(*,*) 'IAM\MUT: iXSN ', iXSN
    !endif

    
    !
    ! update counters
    !
    if (oSlim.eq.0) then

        if (iXSN.ge.0) then

            ! existing change
            if (icn2.eq.0) then
                !write(*,*) 'IAM\MUT: Existing change at this location (2).'
                !write(*,*) 'IAM\MUT: isp1, ia1 = ', isp1, ia1
                !write(*,*) 'IAM\MUT: icn1, int1 = ', icn1, int1
                if (a_cmr(icn1,isp1,ia1,1).eq.1) then
                    a_nm(1,isp1,ia1,g_nt_o(ig1,int1),nt1) = a_nm(1,isp1,ia1,g_nt_o(ig1,int1),nt1) - 1
                    a_nmt(1,isp1,ia1) = a_nmt(1,isp1,ia1) - 1
                    if (a_cXSN(icn1,isp1,ia1,1).eq.1) then
                        a_nmS(1,isp1,ia1) = a_nmS(1,isp1,ia1) - 1
                    elseif (a_cXSN(icn1,isp1,ia1,1).eq.2) then
                        a_nmN(1,isp1,ia1) = a_nmN(1,isp1,ia1) - 1
                    endif
                    !write(*,*) 'IAM\MUT: Undo existing mutation.'
                    !write(*,*) 'IAM\MUT: icn1 = ', icn1
                    !write(*,*) 'IAM\MUT: a_nmt = ', a_nmt(1,isp1,ia1)
                    !write(*,*) 'IAM\MUT: a_nrt = ', a_nrt(1,isp1,ia1)
                elseif (a_cmr(icn1,isp1,ia1,1).eq.2) then
                    a_nr(1,isp1,ia1,g_nt_o(ig1,int1),nt1) = a_nr(1,isp1,ia1,g_nt_o(ig1,int1),nt1) - 1
                    a_nrt(1,isp1,ia1) = a_nrt(1,isp1,ia1) - 1
                    if (a_cXSN(icn1,isp1,ia1,1).eq.1) then
                        a_nrS(1,isp1,ia1) = a_nrS(1,isp1,ia1) - 1
                    elseif (a_cXSN(icn1,isp1,ia1,1).eq.2) then
                        a_nrN(1,isp1,ia1) = a_nrN(1,isp1,ia1) - 1
                    endif
                    !write(*,*) 'IAM\MUT: Undo existing recombination.'
                    !write(*,*) 'IAM\MUT: icn1 = ', icn1
                    !write(*,*) 'IAM\MUT: a_nmt = ', a_nmt(1,isp1,ia1)
                    !write(*,*) 'IAM\MUT: a_nrt = ', a_nrt(1,isp1,ia1) 
                endif
            endif
            
            ! new change
            if (g_nt_o(ig1,int1).ne.nt2) then
                a_nm(1,isp1,ia1,g_nt_o(ig1,int1),nt2) = a_nm(1,isp1,ia1,g_nt_o(ig1,int1),nt2) + 1
                a_nmt(1,isp1,ia1) = a_nmt(1,isp1,ia1) + 1
                if (iXSN.eq.1) then
                    a_nmS(1,isp1,ia1) = a_nmS(1,isp1,ia1) + 1
                elseif (iXSN.eq.2) then
                    a_nmN(1,isp1,ia1) = a_nmN(1,isp1,ia1) + 1
                endif
            endif
            !write(*,*) 'IAM\MUT: Do new mutation.'
            !write(*,*) 'IAM\MUT: icn1 = ', icn1
            !write(*,*) 'IAM\MUT: a_nmt = ', a_nmt(1,isp1,ia1)
            !write(*,*) 'IAM\MUT: a_nrt = ', a_nrt(1,isp1,ia1) 
        endif

        if (iXSN.ge.0) then

            nmTt_sp(isp1) = nmTt_sp(isp1) + 1.
            nmTij_sp(isp1,nt1,nt2) = nmTij_sp(isp1,nt1,nt2) + 1.
            if (iXSN.eq.0) then
                nmXt_sp(isp1) = nmXt_sp(isp1) + 1.
                nmXij_sp(isp1,nt1,nt2) = nmXij_sp(isp1,nt1,nt2) + 1.
            elseif (iXSN.eq.1) then
                nmSt_sp(isp1) = nmSt_sp(isp1) + 1.
                nmSij_sp(isp1,nt1,nt2) = nmSij_sp(isp1,nt1,nt2) + 1.
            else
                nmNt_sp(isp1) = nmNt_sp(isp1) + 1.
                nmNij_sp(isp1,nt1,nt2) = nmNij_sp(isp1,nt1,nt2) + 1.
            endif

        endif
    
    endif

    !if (isp1.eq.1) then
    !if (a_sn(1,isp1,ia1).eq.390) then
    !write(*,*) 'IAM\MUT: done ixmr = ', ixmr
    !endif

    enddo ! ixmr
    
endif ! srm

!
! --- recombination ---
!
if (a_srr(1,isp1,ia1).ge.1.) then

    !if (isp1.eq.4) then
    !write(*,*) 'IAM\REC: t_sp, ia1 = ', t_sp(isp1), ia1
    !endif
    ig1 = a_gid(1,isp1,ia1)
    a_srr(1,isp1,ia1) = 0.

    !if (oHANS.eq.-1) then
    !    xRuse = bucci(dble(xR), dble(xR*10), idum_sp(isp1))
    !else
        xRuse = xR
    !endif

    do ixmr = 1, xRuse
    !if (isp1.eq.4) then
    !write(*,*) 'IAM\REC: ixmr, xRuse = ', ixmr, xRuse
    !endif
    
    !
    ! select donor
    !
    ia4 = 0
401 continue
    if (oHANS.ge.1) then
        r = r4_uni(idum_sp(isp1))
        ia3 = 1 + r * nar_sp(isp1)*xR
        if (ia3.gt.nar_sp(isp1)*xR) then
            ia3 = nar_sp(isp1)*xR
        endif
        ia2 = iad_sp(isp1,ia3)
    else
        r = r4_uni(idum_sp(isp1))
        ia2 = 1 + r * nat_sp(1,isp1)
        if (ia2.gt.nat_sp(1,isp1)) then
            ia2 = nat_sp(1,isp1)
        endif
        if (a_sa(1,isp1,ia2).eq.0) then
            !write(*,*) 'IAM\REC: Inactive agent, trying again...'
            goto 401
        endif
    endif
    if ((ia1.eq.ia2).or.(a_sa(1,isp1,ia2).eq.0)) then
        ia4 = ia4 + 1
        if (ia4.gt.(nar_sp(isp1)*xR)) then
            write(*,*) 'IAM\REC: Warning. Problem finding donor.'
            write(*,*) 'IAM\REC: isp1, nar_sp = ', isp1, nar_sp(isp1)
            goto 403
        endif
        !write(*,*) 'IAM\REC: Selfie, trying again...'
        goto 401
    endif
    ig2 = a_gid(1,isp1,ia2)
    !if (isp1.eq.4) then
    !write(*,*) 'IAM\REC: ia2, ig2 = ', ia2, ig2
    !endif

    !
    ! select length
    !
402 continue
    if (oR.eq.1) then
        lR = nint(lRm)
    else
        r = r4_uni(idum_sp(isp1))
        if (r.eq.0.) then
            goto 402
        else
            lR = nint(-log(r)*lRm)
        endif
    endif
    if ((lR.lt.1).or.(lR.gt.lRmx)) then
        goto 402
    endif
    !if (isp1.eq.4) then
    !write(*,*) 'IAM\REC: lR = ', lR
    !endif
    rdeltat_sp(isp1) = rdeltat_sp(isp1) + dble(lR) * a_sr(1,isp1,ia1)

    !
    ! select start position
    !
    r = r4_uni(idum_sp(isp1))
    Rstart = 1 + r * (dimnnt - lR)
    Rend = Rstart + lR - 1 
    !if (isp1.eq.4) then
    !write(*,*) 'IAM\REC: Rstart, Rend = ', Rstart, Rend
    !endif

    !
    ! update sequences
    !
    !cp_mr2a_t1_sp(isp1) = omp_get_wtime()
    do int1 = max(1,Rstart-2), min(dimnnt,Rend+2)
        g_ntu1(isp1,int1) = g_nt(ig1,int1)
        g_ntu2(isp1,int1) = g_nt(ig2,int1)
        g_cmru1(isp1,int1) = g_cmr(ig1,int1)
        g_cmru2(isp1,int1) = g_cmr(ig2,int1)
        g_cic1(isp1,int1) = -1
        !g_cic2(isp1,int1) = -1
        g_cfp1(isp1,int1) = 1.
        g_cfp2(isp1,int1) = 1.
    enddo
    !cp_mr2a_delta_sp(isp1) = cp_mr2a_delta_sp(isp1) + omp_get_wtime() - cp_mr2a_t1_sp(isp1)
    ! note: profiling suggests these two loops below account for about 80% of CPU time
    !cp_mr2b_t1_sp(isp1) = omp_get_wtime()
    do ic1 = 1, a_cn(1,isp1,ia1)
        int1 = a_ck(ic1,isp1,ia1,1)
        if ((int1.ge.Rstart-2).and.(int1.le.Rend+2)) then
            g_ntu1(isp1,int1) = a_cnt(ic1,isp1,ia1,1)
            g_cic1(isp1,int1) = ic1
            g_cfp1(isp1,int1) = a_cfp(ic1,isp1,ia1,1)
            if (oSlim.eq.0) then
                g_cmru1(isp1,int1) = a_cmr(ic1,isp1,ia1,1)
            endif
        endif
    enddo
    !cp_mr2b_delta_sp(isp1) = cp_mr2b_delta_sp(isp1) + omp_get_wtime() - cp_mr2b_t1_sp(isp1)
    !cp_mr2c_t1_sp(isp1) = omp_get_wtime()
    do ic1 = 1, a_cn(1,isp1,ia2)
        int1 = a_ck(ic1,isp1,ia2,1)
        if ((int1.ge.Rstart-2).and.(int1.le.Rend+2)) then
            g_ntu2(isp1,int1) = a_cnt(ic1,isp1,ia2,1)
            !g_cic2(isp1,int1) = ic1
            g_cfp2(isp1,int1) = a_cfp(ic1,isp1,ia2,1)
            if (oSlim.eq.0) then
                g_cmru2(isp1,int1) = a_cmr(ic1,isp1,ia2,1)
            endif
        endif
    enddo
    !cp_mr2c_delta_sp(isp1) = cp_mr2c_delta_sp(isp1) + omp_get_wtime() - cp_mr2c_t1_sp(isp1)
    
    !
    ! copy
    !
    int4 = 0
    do int1 = Rstart, Rend
        nt1 = g_ntu1(isp1,int1)
        nt2 = g_ntu2(isp1,int1)
        icn1 = g_cic1(isp1,int1)
        fp1 = g_cfp1(isp1,int1)
        fp2 = g_cfp2(isp1,int1)
        if (nt1.ne.nt2) then
            !if (isp1.eq.4) then
            !write(*,*) 'IAM\REC: int1 = ', int1
            !write(*,*) 'IAM\REC: nt1, nt2 = ', nt1, nt2
            !endif
            do int5 = int1 - 2, int1 + 2
                int3 = int5 - int1 + 3
                nt5(int3) = -9
                if ((int5.ge.1).and.(int5.le.dimnnt)) then
                    nt5(int3) = g_ntu1(isp1,int5)
                endif
            enddo
            int4 = int4 + 1
            call rec_c(isp1,ia1,int1,nt5,nt2,2,icn1,iXSN,ierr,fp2,fp1)
            if (ierr.eq.1) then
                write(*,*) 'IAM\REC: ERR 1 in rec_c: jcn > dimnc.', a_cn(1,isp1,ia1), dimnc
                write(*,*) 'IAM\REC: Rstart, Rend = ', Rstart, Rend
                write(*,*) 'IAM\REC: lR, int1 = ', lR, int1
                write(*,*) 'IAM\REC: int4 = ', int4
                write(*,*) 'IAM\REC: isp1, ia1 = ', isp1, ia1
                write(*,*) 'IAM\REC: nt1, nt2 = ', nt1, nt2
                !err_sp(isp1) = 1201
                !goto 9050
                call kill_agent(isp1,ia1)
                ncg_sp(isp1) = ncg_sp_prev
                goto 403
            endif
            if (ierr.eq.2) then
                write(*,*) 'IAM\REC: ERR 2 in rec_c. Invalid iNS.'
                err_sp(isp1) = 1301
                goto 9050
            endif
            if (a_sa(1,isp1,ia1).lt.1) then
                nart_sp(isp1) = nart_sp(isp1) - 1
                nsrrt_sp(isp1) = nsrrt_sp(isp1) - a_sr(1,isp1,ia1)
                nawt_sp(isp1) = nawt_sp(isp1) + 1
                nsrwt_sp(isp1) = nsrwt_sp(isp1) + a_sr(1,isp1,ia1)
                goto 403
            endif

            !g_cic1(isp1,int1) = a_cn(1,isp1,ia1)
            
        endif ! nt1<>nt2

        if (iXSN.lt.0) then
            goto 404
        endif
        
        ! update recombinations count
        if (oSlim.eq.0) then

            !
            ! rtype definition:
            !                   Donor
            !                   None Mut. Reco.
            ! Recipient None    1    2    3
            !           Mut.    4    5    6
            !           Reco.   7    8    9
            !
            !
            if (g_cmru1(isp1,int1).eq.0) then
                if (g_cmru2(isp1,int1).eq.0) then
                    irtype = 1
                elseif (g_cmru2(isp1,int1).eq.1) then
                    irtype = 2
                elseif (g_cmru2(isp1,int1).eq.2) then
                    irtype = 3
                endif
            elseif (g_cmru1(isp1,int1).eq.1) then
                if (g_cmru2(isp1,int1).eq.0) then
                    irtype = 4
                elseif (g_cmru2(isp1,int1).eq.1) then
                    irtype = 5
                elseif (g_cmru2(isp1,int1).eq.2) then
                    irtype = 6
                endif
            elseif (g_cmru1(isp1,int1).eq.2) then
                if (g_cmru2(isp1,int1).eq.0) then
                    irtype = 7
                elseif (g_cmru2(isp1,int1).eq.1) then
                    irtype = 8
                elseif (g_cmru2(isp1,int1).eq.2) then
                    irtype = 9
                endif
            endif

            rtype_sp(isp1,irtype) = rtype_sp(isp1,irtype) + 1

            ! existing mutation
            if (irtype.eq.4) then
                a_nm(1,isp1,ia1,g_nt_o(ig1,int1),nt1) = a_nm(1,isp1,ia1,g_nt_o(ig1,int1),nt1) - 1
                a_nmt(1,isp1,ia1) = a_nmt(1,isp1,ia1) - 1
                if (a_cXSN(icn1,isp1,ia1,1).eq.1) then
                    a_nmS(1,isp1,ia1) = a_nmS(1,isp1,ia1) - 1
                elseif (a_cXSN(icn1,isp1,ia1,1).eq.2) then
                    a_nmN(1,isp1,ia1) = a_nmN(1,isp1,ia1) - 1
                endif
                !write(*,*) '.'
                !write(*,*) 'IAM\REC: Undo existing mutation.'
                !write(*,*) 'IAM\REC: rtype = ', irtype
                !write(*,*) 'IAM\REC: icn1 = ', icn1
                !write(*,*) 'IAM\REC: a_nmt = ', a_nmt(1,isp1,ia1) 
                !write(*,*) 'IAM\REC: a_nrt = ', a_nrt(1,isp1,ia1) 
            endif
            if ((irtype.eq.5).or.(irtype.eq.6)) then
            if (nt1.ne.nt2) then
                a_nm(1,isp1,ia1,g_nt_o(ig1,int1),nt1) = a_nm(1,isp1,ia1,g_nt_o(ig1,int1),nt1) - 1
                a_nmt(1,isp1,ia1) = a_nmt(1,isp1,ia1) - 1
                if (a_cXSN(icn1,isp1,ia1,1).eq.1) then
                    a_nmS(1,isp1,ia1) = a_nmS(1,isp1,ia1) - 1
                elseif (a_cXSN(icn1,isp1,ia1,1).eq.2) then
                    a_nmN(1,isp1,ia1) = a_nmN(1,isp1,ia1) - 1
                endif
                a_nr(1,isp1,ia1,g_nt_o(ig1,int1),nt2) = a_nr(1,isp1,ia1,g_nt_o(ig1,int1),nt2) + 1
                a_nrt(1,isp1,ia1) = a_nrt(1,isp1,ia1) + 1
                if (iXSN.eq.1) then
                    a_nrS(1,isp1,ia1) = a_nrS(1,isp1,ia1) + 1
                elseif (iXSN.eq.2) then
                    a_nrN(1,isp1,ia1) = a_nrN(1,isp1,ia1) + 1
                endif
                !write(*,*) '.'
                !write(*,*) 'IAM\REC: Undo existing mutation.'
                !write(*,*) 'IAM\REC: rtype = ', irtype
                !write(*,*) 'IAM\REC: icn1 = ', icn1
                !write(*,*) 'IAM\REC: a_nmt = ', a_nmt(1,isp1,ia1) 
                !write(*,*) 'IAM\REC: a_nrt = ', a_nrt(1,isp1,ia1) 
            endif
            endif

            ! existing recombination
            if (irtype.eq.7) then
                a_nr(1,isp1,ia1,g_nt_o(ig1,int1),nt1) = a_nr(1,isp1,ia1,g_nt_o(ig1,int1),nt1) - 1
                a_nrt(1,isp1,ia1) = a_nrt(1,isp1,ia1) - 1
                if (a_cXSN(icn1,isp1,ia1,1).eq.1) then
                    a_nrS(1,isp1,ia1) = a_nrS(1,isp1,ia1) - 1
                elseif (a_cXSN(icn1,isp1,ia1,1).eq.2) then
                    a_nrN(1,isp1,ia1) = a_nrN(1,isp1,ia1) - 1
                endif
                !write(*,*) 'IAM\REC: Undo existing recombination.'
                !write(*,*) 'IAM\REC: rtype = ', irtype
                !write(*,*) 'IAM\REC: icn1 = ', icn1
                !write(*,*) 'IAM\REC: a_nmt = ', a_nmt(1,isp1,ia1)
                !write(*,*) 'IAM\REC: a_nrt = ', a_nrt(1,isp1,ia1) 
            endif

            ! new recombination
            if ((irtype.ge.2).and.(irtype.le.3)) then
                if (g_nt_o(ig1,int1).ne.nt2) then
                    a_nr(1,isp1,ia1,g_nt_o(ig1,int1),nt2) = a_nr(1,isp1,ia1,g_nt_o(ig1,int1),nt2) + 1
                    a_nrt(1,isp1,ia1) = a_nrt(1,isp1,ia1) + 1
                    if (iXSN.eq.1) then
                        a_nrS(1,isp1,ia1) = a_nrS(1,isp1,ia1) + 1
                    elseif (iXSN.eq.2) then
                        a_nrN(1,isp1,ia1) = a_nrN(1,isp1,ia1) + 1
                    endif
                endif
                !write(*,*) 'IAM\REC: Do new recombination.'
                !write(*,*) 'IAM\REC: rtype = ', irtype
                !write(*,*) 'IAM\REC: icn1 = ', icn1
                !write(*,*) 'IAM\REC: a_cmr = ', a_cmr(icn1,isp1,ia1,1)
                !write(*,*) 'IAM\REC: a_nmt = ', a_nmt(1,isp1,ia1)
                !write(*,*) 'IAM\REC: a_nrt = ', a_nrt(1,isp1,ia1) 
            endif
            
            if ((nt1.ne.nt2).and.(iXSN.ge.0)) then
                nrTt_sp(isp1) = nrTt_sp(isp1) + 1.
                nrTij_sp(isp1,nt1,nt2) = nrTij_sp(isp1,nt1,nt2) + 1.
                if (iXSN.eq.0) then
                    nrXt_sp(isp1) = nrXt_sp(isp1) + 1.
                    nrXij_sp(isp1,nt1,nt2) = nrXij_sp(isp1,nt1,nt2) + 1.
                elseif (iXSN.eq.1) then
                    nrSt_sp(isp1) = nrSt_sp(isp1) + 1.
                    nrSij_sp(isp1,nt1,nt2) = nrSij_sp(isp1,nt1,nt2) + 1.
                else
                    nrNt_sp(isp1) = nrNt_sp(isp1) + 1.
                    nrNij_sp(isp1,nt1,nt2) = nrNij_sp(isp1,nt1,nt2) + 1.
                endif
            endif
            
        endif
        
404 continue
        
    enddo ! int1
    
    rnut_sp(isp1) = rnut_sp(isp1) + dble(int4) / dble(lR) * a_sr(1,isp1,ia1)
    rn_sp(isp1) = rn_sp(isp1) + a_sr(1,isp1,ia1)

    !if (isp1.eq.4) then
    !write(*,*) 'IAM\REC: done ixmr = ', ixmr
    !endif

    enddo ! ixmr
    
endif ! srr

403 continue
endif ! sa
enddo ! ia1
cp_mr2_delta_sp(isp1) = cp_mr2_delta_sp(isp1) + omp_get_wtime() - cp_mr2_t1_sp(isp1)

VmaxCa_sp(isp1) = VmaxCa_sp(isp1) / SRT_sp(isp1)
VmaxNa_sp(isp1) = VmaxNa_sp(isp1) / SRT_sp(isp1)
VmaxPa_sp(isp1) = VmaxPa_sp(isp1) / SRT_sp(isp1)

!
! ----------------------
! --- fixed pop size ---
! ----------------------
!
if ((oFixCells.eq.1).and.(oHANS.eq.0)) then

PD = 1. - Popic / dble(naa)

if (PD.gt.0.) then
do ia1 = 1, nat_sp(1,isp1)
if (a_sa(1,isp1,ia1).ge.1) then

    r = r4_uni(idum_sp(isp1))
    if (r.le.PD) then
        !write(*,*) 'IAM: Death/outflow. ia1, sn = ', ia1, a_sn
        !nwt_sp(isp1) = nwt_sp(isp1) + 1
        call kill_agent(isp1,ia1)
    endif

endif
enddo
endif

endif

!
! ----------------
! --- dilution ---
! ----------------
!
ipD = 1
if (oHANS.eq.0) then
if (dtD.gt.0.) then
ipD = 0
if ((dtD.gt.0.).and.(t_sp(isp1).ge.tDNext_sp(isp1))) then
    !write(*,*) 'IAM\D: t_sp = ', t_sp(isp1)
    
    ipD = 1
    PD = 1. - fD
    if (oFixCells.eq.2) then
        PS = dble(nic) / dble(nat_sp(1,isp1)-naf_sp(isp1))
        PD = 1. - PS
    endif
    !write(*,*) 'IAM\D: PD = ', PD
    do ia1 = 1, nat_sp(1,isp1)
    if (a_sa(1,isp1,ia1).ge.1) then
        r = r4_uni(idum_sp(isp1))
        if (r.le.PD) then
            !write(*,*) 'IAM\D: t_sp, SN = ', t_sp(isp1), a_sn(1,isp1,ia1)
            call kill_agent(isp1,ia1)
        endif
    endif
    enddo

    SC_sp(isp1) = SC_sp(isp1) * fD + SinC * (1. - fD)
    SN_sp(isp1) = SN_sp(isp1) * fD + SinN * (1. - fD)
    SP_sp(isp1) = SP_sp(isp1) * fD + SinP * (1. - fD)

    tDNext_sp(isp1) = t_sp(isp1) + dtD

endif
endif
endif

!
! -------------------------------
! --- calc extracellular conc ---
! -------------------------------
!
if (oFixExt.eq.0) then

    if (oImpExt.eq.0) then

        SC_sp(isp1) = SC_sp(isp1) + RSC_sp(isp1) * dt
        SN_sp(isp1) = SN_sp(isp1) + RSN_sp(isp1) * dt
        SP_sp(isp1) = SP_sp(isp1) + RSP_sp(isp1) * dt

        if (SC_sp(isp1).lt.0.) then
            write(*,*) 'IAM\ERR: SC < 0. ', t_sp(isp1), SC_sp(isp1)
            err_sp(isp1) = 1401
            goto 9050
        endif
        if (isnan(SC_sp(isp1))) then
            write(*,*) 'IAM\ERR: SC = nan. ', t_sp(isp1), SC_sp(isp1)
            err_sp(isp1) = 1402
            goto 9050
        endif

        if (SN_sp(isp1).lt.0.) then
            write(*,*) 'IAM\ERR: SN < 0. ', t_sp(isp1), SN_sp(isp1)
            err_sp(isp1) = 1403
            goto 9050
        endif
        if (isnan(SN_sp(isp1))) then
            write(*,*) 'IAM\ERR: SN = nan. ', t_sp(isp1), SN_sp(isp1)
            err_sp(isp1) = 1404
            goto 9050
        endif

        if (SP_sp(isp1).lt.0.) then
            write(*,*) 'IAM\ERR: SP < 0. ', t_sp(isp1), SP_sp(isp1)
            err_sp(isp1) = 1405
            goto 9050
        endif
        if (isnan(SP_sp(isp1))) then
            write(*,*) 'IAM\ERR: SP = nan. ', t_sp(isp1), SP_sp(isp1)
            err_sp(isp1) = 1406
            goto 9050
        endif
        
    else
        
        SUMVmax = SRT_sp(isp1) * VmaxCa_sp(isp1)
        SC_sp(isp1) = (Q_sp(isp1)*SinC*dt+VOL_sp(isp1)*SC_sp(isp1)-VOL_sp(isp1)*KmC-Q_sp(isp1)*dt*KmC-SUMVmax*dt+ &
            (Q_sp(isp1)**2*dt**2*KmC**2+2*VOL_sp(isp1)**2*SC_sp(isp1)*KmC+VOL_sp(isp1)**2*SC_sp(isp1)**2+ &
            VOL_sp(isp1)**2*KmC**2+SUMVmax**2*dt**2+2*Q_sp(isp1)*SinC*dt*VOL_sp(isp1)*SC_sp(isp1)+ &
            2*Q_sp(isp1)*SinC*dt*VOL_sp(isp1)*KmC+2*Q_sp(isp1)**2*SinC*dt**2*KmC- &
            2*Q_sp(isp1)*SinC*dt**2*SUMVmax+2*VOL_sp(isp1)*SC_sp(isp1)*Q_sp(isp1)*dt*KmC-2*VOL_sp(isp1)*SC_sp(isp1)*SUMVmax*dt+ &
            2*VOL_sp(isp1)*KmC**2*Q_sp(isp1)*dt+2*VOL_sp(isp1)*KmC*SUMVmax*dt+2*Q_sp(isp1)*dt**2*KmC*SUMVmax+Q_sp(isp1)**2*SinC**2*dt**2)**0.5) &
            / (2*(VOL_sp(isp1)+Q_sp(isp1)*dt))

        SUMVmax = SRT_sp(isp1) * VmaxNa_sp(isp1)
        SN_sp(isp1) = (Q_sp(isp1)*SinN*dt+VOL_sp(isp1)*SN_sp(isp1)-VOL_sp(isp1)*KmN-Q_sp(isp1)*dt*KmN-SUMVmax*dt+ &
            (Q_sp(isp1)**2*dt**2*KmN**2+2*VOL_sp(isp1)**2*SN_sp(isp1)*KmN+VOL_sp(isp1)**2*SN_sp(isp1)**2+ &
            VOL_sp(isp1)**2*KmN**2+SUMVmax**2*dt**2+2*Q_sp(isp1)*SinN*dt*VOL_sp(isp1)*SN_sp(isp1)+ &
            2*Q_sp(isp1)*SinN*dt*VOL_sp(isp1)*KmN+2*Q_sp(isp1)**2*SinN*dt**2*KmN- &
            2*Q_sp(isp1)*SinN*dt**2*SUMVmax+2*VOL_sp(isp1)*SN_sp(isp1)*Q_sp(isp1)*dt*KmN-2*VOL_sp(isp1)*SN_sp(isp1)*SUMVmax*dt+ &
            2*VOL_sp(isp1)*KmN**2*Q_sp(isp1)*dt+2*VOL_sp(isp1)*KmN*SUMVmax*dt+2*Q_sp(isp1)*dt**2*KmN*SUMVmax+Q_sp(isp1)**2*SinN**2*dt**2)**0.5) &
            / (2*(VOL_sp(isp1)+Q_sp(isp1)*dt))

        SUMVmax = SRT_sp(isp1) * VmaxPa_sp(isp1)
        SP_sp(isp1) = (Q_sp(isp1)*SinP*dt+VOL_sp(isp1)*SP_sp(isp1)-VOL_sp(isp1)*KmP-Q_sp(isp1)*dt*KmP-SUMVmax*dt+ &
            (Q_sp(isp1)**2*dt**2*KmP**2+2*VOL_sp(isp1)**2*SP_sp(isp1)*KmP+VOL_sp(isp1)**2*SP_sp(isp1)**2+ &
            VOL_sp(isp1)**2*KmP**2+SUMVmax**2*dt**2+2*Q_sp(isp1)*SinP*dt*VOL_sp(isp1)*SP_sp(isp1)+ &
            2*Q_sp(isp1)*SinP*dt*VOL_sp(isp1)*KmP+2*Q_sp(isp1)**2*SinP*dt**2*KmP- &
            2*Q_sp(isp1)*SinP*dt**2*SUMVmax+2*VOL_sp(isp1)*SP_sp(isp1)*Q_sp(isp1)*dt*KmP-2*VOL_sp(isp1)*SP_sp(isp1)*SUMVmax*dt+ &
            2*VOL_sp(isp1)*KmP**2*Q_sp(isp1)*dt+2*VOL_sp(isp1)*KmP*SUMVmax*dt+2*Q_sp(isp1)*dt**2*KmP*SUMVmax+Q_sp(isp1)**2*SinP**2*dt**2)**0.5) &
            / (2*(VOL_sp(isp1)+Q_sp(isp1)*dt))

    endif

endif

!
! --------------------
! --- advance time ---
! --------------------
!
t_sp(isp1) = t_sp(isp1) + dt
tDNext_sp(isp1) = tDNext_sp(isp1)

!write(*,*) 'IAM: Inner time end: isp1 = ', isp1
enddo ! itx

9050 continue

!write(*,*) 'IAM: Subpopulation end: isp1 = ', isp1
enddo ! isp1
!$OMP END DO
!$OMP END PARALLEL
do isp1 = 1, dimnsp
    if (err_sp(isp1).gt.0) then
        write(*,*) 'IAM:'
        write(*,*) 'IAM: ERR. isp1, err_sp = ', isp1, err_sp(isp1)
        goto 9100
    endif
enddo

cp_partot_delta = cp_partot_delta + omp_get_wtime() - cp_partot_t1

!
! --------------------
! --- advance time ---
! --------------------
!
t = t + dt * dble(ntx)
th = th + dt * dble(ntx)
tDNext = tDNext_sp(1)

!
! ------------------
! --- update gsn ---
! ------------------
!
do isp1 = 1, dimnsp
    if (gsn_sp(isp1).gt.gsn) then
        gsn = gsn_sp(isp1)
        if (gsn.gt.gsnmax) then
            write(*,*) 'IAM\SN: Warning: resetting sn.'
            gsn = 1
        endif
    endif
enddo

!
! --------------------------
! --- mix subpopulations ---
! --------------------------
!
if (oIter.eq.1) then
    write(*,*) 'IAM\ITER: Skipping mixing.'
    goto 610
endif

if (dimnsp.gt.1) then
    write(*,*) 'IAM\MIX: Mixing. t =  ', t

    !
    ! --- check ---
    !
    write(*,*) 'IAM\FIX: Updating fixed changes count (before mix).'
    call fix_changes()

    !
    ! --- 1 > 2 ---
    !
    cp_mix1_t1 = omp_get_wtime()

    ! update SRT
    fkill = 0.
    do isp1 = 1, dimnsp
        !write(*,*) 'IAM\MIX1a: isp1 = ', isp1
        SRT = 0.
        nmx = 10
        do imx1 = 1, nmx
            SRmin(imx1) = 1.e20
            SRmax(imx1) = -1.e20
            iamin(imx1) = -9
            iamax(imx1) = -9
        enddo
        do ia1 = 1, nat_sp(1,isp1)
        if (a_sa(1,isp1,ia1).ge.1) then
            SRT = SRT + a_sr(1,isp1,ia1)
            !write(*,*) 'IAM\MIX1a: ia1, SR = ', ia1, a_sr(1,isp1,ia1)
            do imx1 = 1, nmx
                if (a_sr(1,isp1,ia1).lt.SRmin(imx1)) then
                    do imx2 = nmx-1, imx1, -1
                        SRmin(imx2+1) = SRmin(imx2)
                        iamin(imx2+1) = iamin(imx2)
                    enddo
                    SRmin(imx1)= a_sr(1,isp1,ia1)
                    iamin(imx1)= ia1
                    goto 410
                endif
            enddo
410 continue
            do imx1 = 1, nmx
                if (a_sr(1,isp1,ia1).gt.SRmax(imx1)) then
                    do imx2 = nmx-1, imx1, -1
                        SRmax(imx2+1) = SRmax(imx2)
                        iamax(imx2+1) = iamax(imx2)
                    enddo
                    SRmax(imx1)= a_sr(1,isp1,ia1)
                    iamax(imx1)= ia1
                    goto 420
                endif
            enddo
420 continue
        endif ! sa
        enddo ! ia1
        fkill = fkill + (SRT_sp(isp1)-SRT)/SRT_sp(isp1)
        if (SRT_sp(isp1).lt.SRT) then
            write(*,*) 'IAM\MIX1a: ERR: SRT went up.'
            write(*,*) 'IAM\MIX1a: isp1 = ', isp1
            write(*,*) 'IAM\MIX1a: SRT_sp (old) = ', SRT_sp(isp1)
            write(*,*) 'IAM\MIX1a: SRT (new) = ', SRT
            write(*,*) 'IAM\MIX1a: nat_sp, naf_sp = ', nat_sp(1,isp1), naf_sp(isp1)
            err_sp(isp1) = 1501
            goto 9100
        endif
        SRT_sp(isp1) = SRT
        !write(*,*) 'IAM\MIX1a: isp1, SRT = ', isp1, SRT
        !write(*,*) 'IAM\MIX1a: SRmin, SRmax = ', SRmin, SRmax
        if ((SRmax(1)/SRT).gt.(1./dble(dimnsp))) then
            write(*,*) 'IAM\MIX1a: Warning: SRmax/SRT > 1/dimnsp.'
            write(*,*) 'IAM\MIX1a: isp1, SRT = ', isp1, SRT
            do imx1 = 1, nmx
                write(*,*) 'IAM\MIX1a: imx, ia, SRmax = ', imx1, iamax(imx1), SRmax(imx1)
            enddo
            write(*,*) 'IAM\MIX1a: imx, ia, SRmax = ...'
            do imx1 = nmx, 1, -1
                write(*,*) 'IAM\MIX1a: imx, ia, SRmin = ', imx1, iamin(imx1), SRmin(imx1)
            enddo
        endif
    enddo ! isp1
    write(*,*) 'IAM\MIX1a: fkill = ', fkill/dimnsp
    
    ! handle nat_sp = 0
    nspe = 0
    do isp1 = 1, dimnsp
        !write(*,*) 'IAM\MIX1b: isp1, naa, SRT = ', isp1, (nat_sp(1,isp1)-naf_sp(isp1)), SRT_sp(isp1)
        if ((nat_sp(1,isp1)-naf_sp(isp1)).eq.0) then
            write(*,*) 'IAM\MIX1b: Warning: naa_sp = 0. isp1 = ', isp1
            isp3 = 0
            nspe = nspe + 1
501 continue
            r = r4_uni(idum_sp(isp1))
            isp2 = 1 + r * dimnsp
            if ((nat_sp(1,isp2)-naf_sp(isp2)).eq.0) then
                isp3 = isp3 + 1
                if (isp3.gt.(10*dimnsp)) then
                    write(*,*) 'IAM\MIX1b: ERR. Problem finding good subpopulation.'
                    err_sp(isp1) = 1503
                    goto 9300
                endif
                goto 501
            endif
            write(*,*) 'IAM\MIX1b: Replacing with isp2 = ', isp2
            nat_sp(1,isp1)=nat_sp(1,isp2)
            SRT_sp(isp1)=SRT_sp(isp2)
            do ia1 = 1, nat_sp(1,isp2)
                call copy_agent(1,isp2,ia1,1,isp1,ia1)
            enddo
        endif
    enddo ! isp1
    if (nspe.gt.0) then
        write(*,*) 'IAM\MIX1b: Warning: nat_sp = 0. nspe = ', nspe
    endif
    
    ! spread
    do isp1 = 1, dimnsp
        nat_sp(2,isp1) = 0
    enddo
    do isp1 = 1, dimnsp
        !write(*,*) 'IAM\MIX1c: isp1, SRT_sp = ', isp1, SRT_sp(isp1)
        r = r4_uni(idum_sp(isp1))
        isp2 = 1 + r * dimnsp
        isp3 = 1
        isp4 = 0
        !write(*,*) 'IAM\MIX1c: First isp2, isp3 = ', isp2, isp3
        SRCDF = 0.
        do ia1 = 1, nat_sp(1,isp1)
        if (a_sa(1,isp1,ia1).ge.1) then
            
            if (oMix.eq.1) then
                SRCDF = SRCDF + a_sr(1,isp1,ia1) / SRT_sp(isp1)
                !write(*,*) 'IAM\MIX1c: sr, SRCDF = ', a_sr(1,isp1,ia1), SRCDF
                if (SRCDF.gt.(dble(isp3)/dble(nspmix))) then
                    goto 510
                endif
            else
                if (ia1.gt.(dble(isp3)/dble(nspmix))) then
                    goto 510
                endif
            endif
            goto 520
510 continue
            isp2 = isp2 + 1
            if (isp2.gt.dimnsp) then
                isp2 = 1
            endif
            if (isp3.lt.dimnsp) then
                isp3 = isp3 + 1
            endif
            !write(*,*) 'IAM\MIX1c: isp2, isp3 = ', isp2, isp3
520 continue
            if ((nat_sp(2,isp2)+1).gt.dimna_sp) then
                isp4 = isp4 + 1
                if (isp4.gt.dimnsp) then
                    write(*,*) 'IAM\MIX1c: ERR: all sp at dimna_sp.'
                    err_sp(isp1) = 15051
                    goto 9100
                endif
                goto 510
            endif

            nat_sp(2,isp2)=nat_sp(2,isp2)+1
            if (nat_sp(2,isp2).gt.dimna_sp) then
                write(*,*) 'IAM\MIX1c: ERR: nat_sp > dimna_sp.', nat_sp(2,isp2), dimna_sp
                err_sp(isp1) = 1505
                goto 9100
            endif
            ia2=nat_sp(2,isp2)
            !write(*,*) 'IAM\MIX1c: ia1, ia2 = ', ia1, ia2
            call copy_agent(1,isp1,ia1,2,isp2,ia2)
        endif ! sa>=1
        enddo ! ia1
        !write(*,*) 'IAM\MIX1c: Last isp2, isp3 = ', isp2, isp3

    enddo ! isp1
    cp_mix1_delta = cp_mix1_delta + omp_get_wtime() - cp_mix1_t1
    
    !
    ! --- 2 > 1 ---
    !
    cp_mix2_t1 = omp_get_wtime()

    ! handle nat_sp = 0
    nspe = 0
    do isp1 = 1, dimnsp
        !write(*,*) 'IAM\MIX2a: isp1, nat_sp = ', isp1, nat_sp(2,isp1)
        if (nat_sp(2,isp1).eq.0) then
            !write(*,*) 'IAM\MIX2a: Warning: nat_sp = 0. isp1 = ', isp1
            isp3 = 0
            nspe = nspe + 1
502 continue
            r = r4_uni(idum_sp(isp1))
            isp2 = 1 + r * dimnsp
            if (nat_sp(2,isp2).eq.0) then
                isp3 = isp3 + 1
                if (isp3.gt.(10*dimnsp)) then
                    write(*,*) 'IAM\MIX2a: Problem finding good subpopulation.'
                endif
                goto 502
            endif
            !write(*,*) 'IAM\MIX2a: Replacing with isp2 = ', isp2
            nat_sp(2,isp1)=nat_sp(2,isp2)
            do ia1 = 1, nat_sp(2,isp2)
                call copy_agent(2,isp2,ia1,2,isp1,ia1)
            enddo
        endif
    enddo ! isp1
    if (nspe.gt.0) then
        write(*,*) 'IAM\MIX2a: Warning: nat_sp = 0. nspe = ', nspe
    endif
    
    ! copy
    !$OMP parallel num_threads(ntd) default(shared) PRIVATE(isp1, ia1, ia2)
    !$OMP DO
    do isp1 = 1, dimnsp

        nat_sp(1,isp1) = 0
        naf_sp(isp1) = 0
        SRT_sp(isp1) = 0.
        do ia1 = 1, nat_sp(2,isp1)

            nat_sp(1,isp1)=nat_sp(1,isp1)+1
            if (nat_sp(1,isp1).gt.dimna_sp) then
                write(*,*) 'IAM\MIX2b: ERR: nat_sp > dimna_sp.', nat_sp(1,isp1), dimna_sp
                err_sp(isp1) = 1510
                goto 9060
            endif
            ia2=nat_sp(1,isp1)

            SRT_sp(isp1) = SRT_sp(isp1) + a_sr(2,isp1,ia1)

            call copy_agent(2,isp1,ia1,1,isp1,ia2)
            
        enddo ! ia1

9060 continue
        
    enddo ! isp1
    !$OMP END DO
    !$OMP END PARALLEL
    do isp1 = 1, dimnsp
        if (err_sp(isp1).gt.0) then
            write(*,*) 'IAM:'
            write(*,*) 'IAM: ERR. isp1, err_sp = ', isp1, err_sp(isp1)
            goto 9100
        endif
    enddo
    cp_mix2_delta = cp_mix2_delta + omp_get_wtime() - cp_mix2_t1

endif ! dimnsp>1

SRT = 0.
do isp1 = 1, dimnsp
    SRT = SRT + SRT_sp(isp1)
enddo

610 continue

!do isp1 = 1, dimnsp
!    Popic_sp(isp1) = SRT_sp(isp1) * Popic / SRT
!    write(*,*) 'IAM: isp1, Popic_sp', isp1, Popic_sp(isp1) 
!enddo

!note: these statistics are refering to the pre-mix sp
ncg_sp_mx = 0
nat = 0
nat_sp_mx = 0
naf = 0
SRT = 0.
do isp1 = 1, dimnsp
    if (ncg_sp(isp1).gt.ncg_sp_mx) then
        ncg_sp_mx = ncg_sp(isp1)
    endif
    nmTt = nmTt + nmTt_sp(isp1)
    nmSt = nmSt + nmSt_sp(isp1)
    nmNt = nmNt + nmNt_sp(isp1)
    nmXt = nmXt + nmXt_sp(isp1)
    do nt1 = 1, 4
        do nt2 = 1, 4
            nmTij(nt1,nt2) = nmTij(nt1,nt2) + nmTij_sp(isp1,nt1,nt2)
            nmSij(nt1,nt2) = nmSij(nt1,nt2) + nmSij_sp(isp1,nt1,nt2)
            nmNij(nt1,nt2) = nmNij(nt1,nt2) + nmNij_sp(isp1,nt1,nt2)
            nmXij(nt1,nt2) = nmXij(nt1,nt2) + nmXij_sp(isp1,nt1,nt2)
        enddo
    enddo
    nrTt = nrTt + nrTt_sp(isp1)
    nrSt = nrSt + nrSt_sp(isp1)
    nrNt = nrNt + nrNt_sp(isp1)
    nrXt = nrXt + nrXt_sp(isp1)
    do nt1 = 1, 4
        do nt2 = 1, 4
            nrTij(nt1,nt2) = nrTij(nt1,nt2) + nrTij_sp(isp1,nt1,nt2)
            nrSij(nt1,nt2) = nrSij(nt1,nt2) + nrSij_sp(isp1,nt1,nt2)
            nrNij(nt1,nt2) = nrNij(nt1,nt2) + nrNij_sp(isp1,nt1,nt2)
            nrXij(nt1,nt2) = nrXij(nt1,nt2) + nrXij_sp(isp1,nt1,nt2)
        enddo
    enddo
    rdeltat = rdeltat + rdeltat_sp(isp1)
    rnut = rnut + rnut_sp(isp1)
    rn = rn + rn_sp(isp1)
    do int1 = 1, 9
        rtype(int1) = rtype(int1) + rtype_sp(isp1, int1)
        rtypet = rtypet + rtype_sp(isp1, int1)
    enddo
    fpt = fpt + fpt_sp(isp1)
    fpn2 = fpn2 + fpn2_sp(isp1)
    nat = nat + nat_sp(1,isp1)
    if (nat_sp(1,isp1).gt.nat_sp_mx) then
        nat_sp_mx = nat_sp(1,isp1)
    endif
    naf = naf + naf_sp(isp1)
    nadt = nadt + nadt_sp(isp1)
    nsrdt = nsrdt + nsrdt_sp(isp1)
    fmutt = fmutt + fmutt_sp(isp1)
    frect = frect + frect_sp(isp1)
    fmutrecn = fmutrecn + fmutrecn_sp(isp1)
    namt = namt + namt_sp(isp1)
    nsrmt = nsrmt + nsrmt_sp(isp1)
    nart = nart + nart_sp(isp1)
    nsrrt = nsrrt + nsrrt_sp(isp1)
    nawt = nawt + nawt_sp(isp1)
    nsrwt = nsrwt + nsrwt_sp(isp1)
    SRT = SRT + SRT_sp(isp1)
    if (err_sp(isp1).gt.0) then
        write(*,*) 'IAM:'
        write(*,*) 'IAM: ERR. isp1, err_sp = ', isp1, err_sp(isp1)
        goto 9100
    endif
enddo

!
! ---------------------
! --- check for end ---
! ---------------------
!
naa = nat - naf
if (naa.eq.0) then
    write(*,*) 'IAM:'
    write(*,*) 'IAM: Exit on no agents. t = ', t
    goto 9100
endif

if (t.gt.tend) then
    write(*,*) 'IAM:'
    write(*,*) 'IAM: Normal end. t, tend = ', t, tend
    !goto 9100
    ipSTAf = 1
    goto 5210
endif

!write(*,*) 'IAM: Outer time end.'
enddo ! outer time loop

!
! -------------
! --- CHECK ---
! -------------
!
!do isp1 = 1, dimnsp
!do ia1 = 1, nat_sp(1,isp1)
!    cmt = 0
!    crt = 0
!    if (a_sa(1,isp1,ia1).ge.1) then
!        do ic1 = 1, a_cn(1,isp1,ia1)
!            if (a_cmr(ic1,isp1,ia1,1).eq.1) then
!                cmt = cmt + 1
!            elseif (a_cmr(ic1,isp1,ia1,1).eq.2) then
!                crt = crt + 1
!            endif
!        enddo
!    endif
!    if (cmt.ne.a_nmt(1,isp1,ia1)) then
!        write(*,*) 'IAM\ERR1: ', isp1, ia1
!        write(*,*) 'IAM\ERR1: ', a_sn(1,isp1,ia1)
!        write(*,*) 'IAM\ERR1: ', cmt, a_nmt(1,isp1,ia1)
!        goto 9100
!    endif
!    if (crt.ne.a_nrt(1,isp1,ia1)) then
!        write(*,*) 'IAM\ERR2: ', isp1, ia1
!        write(*,*) 'IAM\ERR2: ', a_sn(1,isp1,ia1)
!        write(*,*) 'IAM\ERR2: ', crt, a_nrt(1,isp1,ia1)
!        goto 9100
!    endif
!enddo
!enddo

9100 continue


!
! ---------------------------
! --- write hotstart file ---
! ---------------------------
!
if ((oHOT.eq.1).or.(oHOT.eq.3)) then
    write(*,*) 'IAM\HOT: '
    write(*,*) 'IAM\HOT: Writing hotstart files.'

    !
    ! --- population ---
    !
    ! mark for check
    ! note: srm and srr variables are used temporarily to track a few cells from model write (here)
    ! to mix read to mix write to model read. srm = flag, srr icomp (assigned in mix)
    npCHK = 10
    do ia1 = 1, npCHK
9120 continue
        r = r4_uni(jdum)
        isp1 = 1 + r * dimnsp
        r = r4_uni(jdum)
        ia2 = 1 + r * nat_sp(1,isp1)
        if(a_sa(1,isp1,ia2).eq.0) then
            goto 9120
        endif 
        a_srm(1,isp1,ia2) = -9.
        write(*,*) 'IAM\HOT\CHK: out sn = ', a_sn(1,isp1,ia2)
    enddo
    !write
    if (cTime.gt.0) then
        open(421,file='IAM_R_POP'//cString//'.bin', form='unformatted')
    else
        open(421,file='IAM_R_POP.bin', form='unformatted')
    endif
    natHOT = 0
    SRTHOT = 0.
    do isp1 = 1, dimnsp
    ia2 = 0
    do ia1 = 1, nat_sp(1,isp1)
    if (a_sa(1,isp1,ia1).ge.1) then
        natHOT = natHOT + 1
        SRTHOT = SRTHOT + a_sr(1,isp1,ia1)
        ia2 = ia2 + 1
        a_srr(1,isp1,ia1) = -9.
        write(421) &
            isp1, &
            ia2, & !< different from read
            a_sa(1,isp1,ia1), &
            a_sr(1,isp1,ia1), &
            a_srm(1,isp1,ia1), & 
            a_srr(1,isp1,ia1), &
            a_iar(1,isp1,ia1), &
            a_gid(1,isp1,ia1), &
            a_cn(1,isp1,ia1)
        write(421) (a_ck(ic1,isp1,ia1,1), ic1 = 1, a_cn(1,isp1,ia1))
        write(421) (a_cnt(ic1,isp1,ia1,1), ic1 = 1, a_cn(1,isp1,ia1))
        do ic1 = 1, a_cn(1,isp1,ia1)
            if ((a_cnt(ic1,isp1,ia1,1).lt.1).or.(a_cnt(ic1,isp1,ia1,1).gt.4)) then
                write(*,*) 'IAM\HOT: ERR. Invalid nt write.'
                write(*,*) 'IAM\HOT: isp1, ia1, ic1 = ', isp1, ia1, ic1
                write(*,*) 'IAM\HOT: nt = ', a_cnt(ic1,isp1,ia1,1)
                goto 9200
            endif
        enddo
        write(421) (a_ntnk(1,isp1,ia1,nt1), nt1 = 1, 4)
        write(421) (a_naa(1,isp1,ia1,iaa1), iaa1 = 1, 20)
        do iaa1 = 1, 20
            write(421) (a_naax(1,isp1,ia1,iaa1,iaa2), iaa2 = 1, 20)
        enddo
        write(421) &
            a_uEg(1, isp1,ia1), &
            a_q0DNAC(1,isp1,ia1), &
            a_q0DNAN(1,isp1,ia1), &
            a_q0DNAP(1,isp1,ia1), &
            a_q0aaC(1,isp1,ia1), &
            a_q0aaN(1,isp1,ia1), &
            a_VmaxC(1,isp1,ia1), &
            a_VmaxN(1,isp1,ia1), &
            a_VmaxP(1,isp1,ia1)
        if (oSlim.eq.0) then
            write(421) &
                a_sn(1,isp1,ia1), &
                a_rsn(1,isp1,ia1), &
                a_nd(1,isp1,ia1), &
                (a_ntnSk(1,isp1,ia1,nt1), nt1 = 1, 4), &
                (a_ntnNk(1,isp1,ia1,nt1), nt1 = 1, 4), &
                (a_ntnXk(1,isp1,ia1,nt1), nt1 = 1, 4), &
                a_nmt(1,isp1,ia1), &
                a_nmS(1,isp1,ia1), &
                a_nmN(1,isp1,ia1), &
                a_nrt(1,isp1,ia1), &
                a_nrS(1,isp1,ia1), &
                a_nrN(1,isp1,ia1)
            write(421) (a_cXSN(ic1,isp1,ia1,1), ic1 = 1, a_cn(1,isp1,ia1))
            write(421) (a_cmr(ic1,isp1,ia1,1), ic1 = 1, a_cn(1,isp1,ia1))
            if ((a_cXSN(ic1,isp1,ia1,1).lt.0).or.(a_cmr(ic1,isp1,ia1,1).lt.0)) then
                write(*,*) 'IAM\HOT: ERR. Invalid cXSN or cmr write.'
                write(*,*) 'IAM\HOT: isp1, ia1, ic1 = ', isp1, ia1, ic1
                write(*,*) 'IAM\HOT: cXSN = ', a_cXSN(ic1,isp1,ia1,1)
                write(*,*) 'IAM\HOT: cmr = ', a_cmr(ic1,isp1,ia1,1)
                goto 9200
            endif
            write(421) ((a_nm(1,isp1,ia1,nt1,nt2), nt1 = 1, 4), nt2 = 1, 4)
            write(421) ((a_nr(1,isp1,ia1,nt1,nt2), nt1 = 1, 4), nt2 = 1, 4)
            write(421) (a_nntt(1,isp1,ia1,nt1), nt1 = 1, 4)
            write(421)  &
                a_VC(1,isp1,ia1), &
                a_VN(1,isp1,ia1), &
                a_VP(1,isp1,ia1), &
                a_qC(1,isp1,ia1), &
                a_qN(1,isp1,ia1), &
                a_qP(1,isp1,ia1), &
                a_CLimit(1,isp1,ia1), &
                a_NLimit(1,isp1,ia1), &
                a_PLimit(1,isp1,ia1), &
                a_tb(1,isp1,ia1), &
                a_tg(1,isp1,ia1)
        endif
    endif ! sa
    enddo ! ia1
    write(*,*) 'IAM\HOT: isp1, ia2 = ', isp1, ia2
    enddo ! isp1
    close(421)
    write(*,*) 'IAM\HOT: nat, SRT = ', natHOT, SRTHOT

    !
    ! --- other ---
    !
    tHot = th
    if (cTime.gt.0) then
        open(411,file='IAM_R_OTH'//cString//'.txt')
    else
        open(411,file='IAM_R_OTH.txt')
    endif
    write(411,*) iHOT, tHOT, natHOT, SRTHOT
    close(411)
    
    !
    ! --- genome ---
    !
    ig1 = 1
    if (ngt.gt.1) then
        write(*,*) 'IAM\HOT: Warning: Multiple genomes. ngt = ', ngt
    endif
    if (g_ncff(ig1).gt.0) then
        write(*,*) 'IAM\HOT: Warning: Flushed changes. g_ncff = ', g_ncff(ig1)
    endif
    if (cTime.gt.0) then
        open(431,file='IAM_R_DNA'//cString//'.txt')
    else
        open(431,file='IAM_R_DNA.txt')
    endif
    write(431,205) '>IAMHOT'//sap(0)//sap(0)//' IAM GENOME HOT ip: '//sap(0)//' ia: '//sap(0)//'.'
    write(431,201) (int_snt(g_nt(1,int1)),int1=1,dimnnt)
    close(431)
    
endif

9200 continue

!
! --- iteration ---
!
if (oIter.eq.1) then
    write(*,*) 'IAM\ITER: Done iIter = ', iIter
    iIter = iIter + 1
    if (iIter.le.nIter) then
        goto 8100
    endif
endif

     
!
! ------------------------------
! --- write final statistics ---
! ------------------------------
!
! --- dimensions ---
!
write(*,*) 'IAM:'
write(*,*) 'IAM: Dimensions:'
write(*,*) 'IAM: dimna_sp, nat_sp_mx = ', dimna_sp, nat_sp_mx
write(*,*) 'IAM: dimng, ngt = ', dimng, ngt
write(*,*) 'IAM: dimnc, ncg_sp_mx = ', dimnc, ncg_sp_mx

!
! --- divisions etc counters ---
!
write(*,*) 'IAM:'
write(*,*) 'IAM: Division etc counters:'
write(*,*) 'IAM: nadt, nsrdt = ', nadt, nsrdt
write(*,*) 'IAM: fmut, frec = ', fmutt/fmutrecn, frect/fmutrecn
write(*,*) 'IAM: nawt, nsrwt = ', nawt, nsrwt
write(*,*) 'IAM: namt, nsrmt = ', namt, nsrmt
write(*,*) 'IAM: nart, nsrrt = ', nart, nsrrt
write(*,*) 'IAM: rdelta, rnu = ', rdelta, rnu
int1 = 1
rtypet = rtypet - rtype(int1)
if (rtypet.gt.0.) then
write(*,*) 'IAM: rtype (i, n, fraction):'
write(*,*) 'IAM:', int1, rtype(int1), 'n/a'
do int1 = 2, 9
write(*,*) 'IAM:', int1, rtype(int1), rtype(int1)/rtypet
enddo
endif
write(*,*) 'IAM: fp = ', fp

write(*,*) 'IAM:'    
write(*,*) 'IAM: Mutation instances'
write(*,*) 'IAM: nmTt n, f = ', nmTt, nmTt/nmTt
write(*,*) 'IAM: nmSt n, f = ', nmSt, nmSt/nmTt
write(*,*) 'IAM: nmNt n, f = ', nmNt, nmNt/nmTt
write(*,*) 'IAM: nmXt n, f = ', nmXt, nmXt/nmTt
do nt1 = 1, 4
    do nt2 = 1,4
        if (nt1.ne.nt2) then
            write(*,*) 'IAM: nt1, nt2  = ', nt1, nt2
            write(*,*) 'IAM: nmTij n, f = ', nmTij(nt1,nt2), nmTij(nt1,nt2)/nmTij(nt1,nt2)
            write(*,*) 'IAM: nmSij n, f = ', nmSij(nt1,nt2), nmSij(nt1,nt2)/nmTij(nt1,nt2)
            write(*,*) 'IAM: nmNij n, f = ', nmNij(nt1,nt2), nmNij(nt1,nt2)/nmTij(nt1,nt2)
            write(*,*) 'IAM: nmXij n, f = ', nmXij(nt1,nt2), nmXij(nt1,nt2)/nmTij(nt1,nt2)
        endif
    enddo
enddo
            
write(*,*) 'IAM\GEN: AT>GC'
int1 = 1
nt11 = 1
nt21 = 4
nt12 = 2
nt22 = 3
bS(int1) = nmSij(nt11,nt21)+nmSij(nt12,nt22)
bN(int1) = nmNij(nt11,nt21)+nmNij(nt12,nt22)
bX(int1) = nmXij(nt11,nt21)+nmXij(nt12,nt22)
bSX(int1) = bS(int1) + bX(int1)
bT(int1) = bN(int1) + bSX(int1)
write(*,*) 'IAM: bS  = ', bS(int1)
write(*,*) 'IAM: bN  = ', bN(int1)
write(*,*) 'IAM: bX  = ', bX(int1)
write(*,*) 'IAM: bSX = ', bSX(int1)
write(*,*) 'IAM: bT  = ', bT(int1)
write(*,*) 'IAM\GEN: GC>AT'
int1 = 2
nt11 = 4
nt21 = 1
nt12 = 3
nt22 = 2
bS(int1) = nmSij(nt11,nt21)+nmSij(nt12,nt22)
bN(int1) = nmNij(nt11,nt21)+nmNij(nt12,nt22)
bX(int1) = nmXij(nt11,nt21)+nmXij(nt12,nt22)
bSX(int1) = bS(int1) + bX(int1)
bT(int1) = bN(int1) + bSX(int1)
write(*,*) 'IAM: bS  = ', bS(int1)
write(*,*) 'IAM: bN  = ', bN(int1)
write(*,*) 'IAM: bX  = ', bX(int1)
write(*,*) 'IAM: bSX = ', bSX(int1)
write(*,*) 'IAM: bT  = ', bT(int1)
write(*,*) 'IAM\GEN: AT>TA'
int1 = 3
nt11 = 1
nt21 = 2
nt12 = 2
nt22 = 1
bS(int1) = nmSij(nt11,nt21)+nmSij(nt12,nt22)
bN(int1) = nmNij(nt11,nt21)+nmNij(nt12,nt22)
bX(int1) = nmXij(nt11,nt21)+nmXij(nt12,nt22)
bSX(int1) = bS(int1) + bX(int1)
bT(int1) = bN(int1) + bSX(int1)
write(*,*) 'IAM: bS  = ', bS(int1)
write(*,*) 'IAM: bN  = ', bN(int1)
write(*,*) 'IAM: bX  = ', bX(int1)
write(*,*) 'IAM: bSX = ', bSX(int1)
write(*,*) 'IAM: bT  = ', bT(int1)
write(*,*) 'IAM\GEN: GC>TA'
int1 = 4
nt11 = 4
nt21 = 2
nt12 = 3
nt22 = 1
bS(int1) = nmSij(nt11,nt21)+nmSij(nt12,nt22)
bN(int1) = nmNij(nt11,nt21)+nmNij(nt12,nt22)
bX(int1) = nmXij(nt11,nt21)+nmXij(nt12,nt22)
bSX(int1) = bS(int1) + bX(int1)
bT(int1) = bN(int1) + bSX(int1)
write(*,*) 'IAM: bS  = ', bS(int1)
write(*,*) 'IAM: bN  = ', bN(int1)
write(*,*) 'IAM: bX  = ', bX(int1)
write(*,*) 'IAM: bSX = ', bSX(int1)
write(*,*) 'IAM: bT  = ', bT(int1)
write(*,*) 'IAM\GEN: AT>CG'
int1 = 5
nt11 = 1
nt21 = 3
nt12 = 2
nt22 = 4
bS(int1) = nmSij(nt11,nt21)+nmSij(nt12,nt22)
bN(int1) = nmNij(nt11,nt21)+nmNij(nt12,nt22)
bX(int1) = nmXij(nt11,nt21)+nmXij(nt12,nt22)
bSX(int1) = bS(int1) + bX(int1)
bT(int1) = bN(int1) + bSX(int1)
write(*,*) 'IAM: bS  = ', bS(int1)
write(*,*) 'IAM: bN  = ', bN(int1)
write(*,*) 'IAM: bX  = ', bX(int1)
write(*,*) 'IAM: bSX = ', bSX(int1)
write(*,*) 'IAM: bT  = ', bT(int1)
write(*,*) 'IAM\GEN: GC>CG'
int1 = 6
nt11 = 4
nt21 = 3
nt12 = 3
nt22 = 4
bS(int1) = nmSij(nt11,nt21)+nmSij(nt12,nt22)
bN(int1) = nmNij(nt11,nt21)+nmNij(nt12,nt22)
bX(int1) = nmXij(nt11,nt21)+nmXij(nt12,nt22)
bSX(int1) = bS(int1) + bX(int1)
bT(int1) = bN(int1) + bSX(int1)
write(*,*) 'IAM: bS  = ', bS(int1)
write(*,*) 'IAM: bN  = ', bN(int1)
write(*,*) 'IAM: bX  = ', bX(int1)
write(*,*) 'IAM: bSX = ', bSX(int1)
write(*,*) 'IAM: bT  = ', bT(int1)
write(*,*) 'IAM\GEN: A+T>G+C'
int1 = 1
int2 = 5
write(*,*) 'IAM\GEN: bS  = ', bS(int1)+bS(int2)
write(*,*) 'IAM\GEN: bN  = ', bN(int1)+bN(int2)
write(*,*) 'IAM\GEN: bX  = ', bX(int1)+bX(int2)
write(*,*) 'IAM\GEN: bSX = ', bSX(int1)+bSX(int2)
write(*,*) 'IAM\GEN: bT  = ', bT(int1)+bT(int2)
write(*,*) 'IAM\GEN: G+C>A+T'
int1 = 2
int2 = 4
write(*,*) 'IAM\GEN: bS  = ', bS(int1)+bS(int2)
write(*,*) 'IAM\GEN: bN  = ', bN(int1)+bN(int2)
write(*,*) 'IAM\GEN: bX  = ', bX(int1)+bX(int2)
write(*,*) 'IAM\GEN: bSX = ', bSX(int1)+bSX(int2)
write(*,*) 'IAM\GEN: bT  = ', bT(int1)+bT(int2)

!
! --- change statistics ---
!
write(*,*) 'IAM:'
write(*,*) 'IAM: Change statisics. Slope is fraction/year.'
write(*,*) 'IAM: t first, last        = ', t_first, t_last
write(*,*) 'IAM: GC first, last       = ', GC_first, GC_last
write(*,*) 'IAM: GCT first, last      = ', GCT_first, GCT_last
write(*,*) 'IAM: GCS first, last      = ', GCS_first, GCS_last
write(*,*) 'IAM: GCN first, last      = ', GCN_first, GCN_last
write(*,*) 'IAM: GCX first, last      = ', GCX_first, GCX_last
write(*,*) 'IAM: GCXS first, last     = ', GCXS_first, GCXS_last
write(*,*) 'IAM: piS(low) first, last = ', piS_first, piS_last
write(*,*) 'IAM: piN(low) first, last = ', piN_first, piN_last
write(*,*) 'IAM: piM(low) first, last = ', piM_first, piM_last
write(*,*) 'IAM: piR(low) first, last = ', piR_first, piR_last
write(*,*) 'IAM: piT(low) first, last = ', piT_first, piT_last
write(*,*) 'IAM: dNdS first, last     = ', dNdS_first, dNdS_last
write(*,*) 'IAM: rm first, last       = ', rm_first, rm_last
write(*,*) 'IAM: q0DNAC first, last   = ', q0DNAC_first, q0DNAC_last
write(*,*) 'IAM: q0DNAN first, last   = ', q0DNAN_first, q0DNAN_last
write(*,*) 'IAM: q0DNAP first, last   = ', q0DNAP_first, q0DNAP_last
write(*,*) 'IAM: q0aaC first, last    = ', q0aaC_first, q0aaC_last
write(*,*) 'IAM: q0aaN first, last    = ', q0aaN_first, q0aaN_last
write(*,*) 'IAM: q0C first, last      = ', q0C_first, q0C_last
write(*,*) 'IAM: q0N first, last      = ', q0N_first, q0N_last
write(*,*) 'IAM: q0P first, last      = ', q0P_first, q0P_last
write(*,*) 'IAM: kg first, last       = ', kg_first, kg_last
if (t_first.eq.t_last) then
    write(*,*) 'IAM: Only one value, skipping slopes.'
else
    write(*,*) 'IAM: GC slope       = ', (GC_last-GC_first)/(t_last-t_first)*365.25
    write(*,*) 'IAM: GCT slope      = ', (GCT_last-GCT_first)/(t_last-t_first)*365.25
    write(*,*) 'IAM: GCS slope      = ', (GCS_last-GCS_first)/(t_last-t_first)*365.25
    write(*,*) 'IAM: GCN slope      = ', (GCN_last-GCN_first)/(t_last-t_first)*365.25
    write(*,*) 'IAM: GCX slope      = ', (GCX_last-GCX_first)/(t_last-t_first)*365.25
    write(*,*) 'IAM: GCXS slope     = ', (GCXS_last-GCXS_first)/(t_last-t_first)*365.25
    write(*,*) 'IAM: piS(low) slope = ', (piS_last-piS_first)/(t_last-t_first)*365.25
    write(*,*) 'IAM: piN(low) slope = ', (piN_last-piN_first)/(t_last-t_first)*365.25
    write(*,*) 'IAM: piM(low) slope = ', (piM_last-piM_first)/(t_last-t_first)*365.25
    write(*,*) 'IAM: piR(low) slope = ', (piR_last-piR_first)/(t_last-t_first)*365.25
    write(*,*) 'IAM: piT(low) slope = ', (piT_last-piT_first)/(t_last-t_first)*365.25
    write(*,*) 'IAM: dNdS slope     = ', (dNdS_last-dNdS_first)/(t_last-t_first)*365.25
    write(*,*) 'IAM: rm slope       = ', (rm_last-rm_first)/(t_last-t_first)*365.25
    write(*,*) 'IAM: q0DNAC slope   = ', (q0DNAC_last-q0DNAC_first)/(t_last-t_first)*365.25
    write(*,*) 'IAM: q0DNAN slope   = ', (q0DNAN_last-q0DNAN_first)/(t_last-t_first)*365.25
    write(*,*) 'IAM: q0DNAP slope   = ', (q0DNAP_last-q0DNAP_first)/(t_last-t_first)*365.25
    write(*,*) 'IAM: q0aaC slope    = ', (q0aaC_last-q0aaC_first)/(t_last-t_first)*365.25
    write(*,*) 'IAM: q0aaN slope    = ', (q0aaN_last-q0aaN_first)/(t_last-t_first)*365.25
    write(*,*) 'IAM: q0C slope      = ', (q0C_last-q0C_first)/(t_last-t_first)*365.25
    write(*,*) 'IAM: q0N slope      = ', (q0N_last-q0N_first)/(t_last-t_first)*365.25
    write(*,*) 'IAM: q0P slope      = ', (q0P_last-q0P_first)/(t_last-t_first)*365.25
    write(*,*) 'IAM: kg slope       = ', (kg_last-kg_first)/(t_last-t_first)*365.25
endif


!
! --------------------------
! --- close output files ---
! --------------------------
!
close(241)
close(242)
close(281)
close(251)

!
! --- code performance ---
!
cp_delta = cp_delta + omp_get_wtime() - cp_t1
call cp
    
!
! ---------------
! ---------------
! --- end msg ---
! ---------------
! ---------------
!
9300 continue
     
write(*,*) 'IAM:'
write(*,*) 'IAM: IAM END'
write(*,*) 'IAM:'

!
! -------------------
! -------------------
! --- subroutines ---
! -------------------
! -------------------
!
contains

!
! ------------------------
! --- code performance ---
! ------------------------
!
subroutine cp()
    cp_ENS_delta = 0.
    cp_HANS1_delta = 0.
    cp_HANS2_delta = 0.
    !cp_HANS2a_delta = 0.
    !cp_HANS2b_delta = 0.
    !cp_HANS2c_delta = 0.
    cp_mr1_delta = 0.
    cp_mr2_delta = 0.
    !cp_mr2a_delta = 0.
    !cp_mr2b_delta = 0.
    !cp_mr2c_delta = 0.
    do isp1 = 1, dimnsp
        cp_ENS_delta = cp_ENS_delta + cp_ENS_delta_sp(isp1)
        cp_HANS1_delta = cp_HANS1_delta + cp_HANS1_delta_sp(isp1)
        cp_HANS2_delta = cp_HANS2_delta + cp_HANS2_delta_sp(isp1)
        !cp_HANS2a_delta = cp_HANS2a_delta + cp_HANS2a_delta_sp(isp1)
        !cp_HANS2b_delta = cp_HANS2b_delta + cp_HANS2b_delta_sp(isp1)
        !cp_HANS2c_delta = cp_HANS2c_delta + cp_HANS2c_delta_sp(isp1)
        cp_mr1_delta = cp_mr1_delta + cp_mr1_delta_sp(isp1)
        cp_mr2_delta = cp_mr2_delta + cp_mr2_delta_sp(isp1)
        !cp_mr2a_delta = cp_mr2a_delta + cp_mr2a_delta_sp(isp1)
        !cp_mr2b_delta = cp_mr2b_delta + cp_mr2b_delta_sp(isp1)
        !cp_mr2c_delta = cp_mr2c_delta + cp_mr2c_delta_sp(isp1)
        !write(*,*) 'IAM\CP: isp1:   ', isp1
        !write(*,*) 'IAM\CP: ENS:    ', cp_ENS_delta_sp(isp1)
        !write(*,*) 'IAM\CP: HANS1:  ', cp_HANS1_delta_sp(isp1)
        !write(*,*) 'IAM\CP: HANS2:  ', cp_HANS2_delta_sp(isp1)
        !write(*,*) 'IAM\CP: HANS2a: ', cp_HANS2a_delta_sp(isp1)
        !write(*,*) 'IAM\CP: HANS2b: ', cp_HANS2b_delta_sp(isp1)
        !write(*,*) 'IAM\CP: HANS2c: ', cp_HANS2c_delta_sp(isp1)
        !write(*,*) 'IAM\CP: mr1:    ', cp_mr1_delta_sp(isp1)
        !write(*,*) 'IAM\CP: mr2:    ', cp_mr2_delta_sp(isp1)
        !write(*,*) 'IAM\CP: mr2a:   ', cp_mr2a_delta_sp(isp1)
        !write(*,*) 'IAM\CP: mr2b:   ', cp_mr2b_delta_sp(isp1)
        !write(*,*) 'IAM\CP: mr2c:   ', cp_mr2c_delta_sp(isp1)
    enddo
    cp_ENS_delta = cp_ENS_delta / dble(dimnsp)
    cp_HANS1_delta = cp_HANS1_delta / dble(dimnsp)
    cp_HANS2_delta = cp_HANS2_delta / dble(dimnsp)
    !cp_HANS2a_delta = cp_HANS2a_delta / dble(dimnsp)
    !cp_HANS2b_delta = cp_HANS2b_delta / dble(dimnsp)
    !cp_HANS2c_delta = cp_HANS2c_delta / dble(dimnsp)
    cp_mr1_delta = cp_mr1_delta / dble(dimnsp)
    cp_mr2_delta = cp_mr2_delta / dble(dimnsp)
    !cp_mr2a_delta = cp_mr2a_delta / dble(dimnsp)
    !cp_mr2b_delta = cp_mr2b_delta / dble(dimnsp)
    !cp_mr2c_delta = cp_mr2c_delta / dble(dimnsp)
    write(*,*) 'IAM\CP:'
    write(*,*) 'IAM\CP: Code performance. Time(h), fraction.'
    write(*,*) 'IAM\CP:'
    write(*,10) 'IAM\CP: Main Loop: ', cp_delta/3600., ', ', cp_delta/cp_delta
    write(*,10) 'IAM\CP: init:      ', cp_init_delta/3600., ', ', cp_init_delta/cp_delta
    write(*,10) 'IAM\CP: iter:      ', cp_iter_delta/3600., ', ', cp_iter_delta/cp_delta
    write(*,10) 'IAM\CP: fix:       ', cp_fix_delta/3600., ', ', cp_fix_delta/cp_delta
    write(*,10) 'IAM\CP: genomes:   ', cp_genomes_delta/3600., ', ', cp_genomes_delta/cp_delta
    write(*,10) 'IAM\CP: output:    ', cp_output_delta/3600., ', ', cp_output_delta/cp_delta
    write(*,10) 'IAM\CP: Par Loop:  ', cp_partot_delta/3600., ', ', cp_partot_delta/cp_delta
    write(*,10) 'IAM\CP: mix1:      ', cp_mix1_delta/3600., ', ', cp_mix1_delta/cp_delta
    write(*,10) 'IAM\CP: mix2(p):   ', cp_mix2_delta/3600., ', ', cp_mix2_delta/cp_delta
    cp_other_delta = cp_delta - ( cp_init_delta + cp_iter_delta + cp_fix_delta + cp_genomes_delta + cp_output_delta + cp_partot_delta + &
                            + cp_mix1_delta + cp_mix2_delta )
    write(*,10) 'IAM\CP: Other:     ', cp_other_delta/3600., ', ', cp_other_delta/cp_delta

    write(*,*) 'IAM\CP:'
    write(*,10) 'IAM\CP: Par Loop:  ', cp_partot_delta/3600., ', ', cp_partot_delta/cp_partot_delta
    write(*,10) 'IAM\CP: ENS:       ', cp_ENS_delta/3600., ', ', cp_ENS_delta/cp_partot_delta
    write(*,10) 'IAM\CP: HANS1:     ', cp_HANS1_delta/3600., ', ', cp_HANS1_delta/cp_partot_delta
    write(*,10) 'IAM\CP: HANS2:     ', cp_HANS2_delta/3600., ', ', cp_HANS2_delta/cp_partot_delta
    write(*,10) 'IAM\CP: mr1:       ', cp_mr1_delta/3600., ', ', cp_mr1_delta/cp_partot_delta
    write(*,10) 'IAM\CP: mr2:       ', cp_mr2_delta/3600., ', ', cp_mr2_delta/cp_partot_delta
    cp_parwait_delta = cp_partot_delta - (cp_ENS_delta + &
                       cp_HANS1_delta + cp_HANS2_delta + &
                       cp_mr1_delta + cp_mr2_delta)
    write(*,10) 'IAM\CP: Wait:      ', cp_parwait_delta/3600., ', ', cp_parwait_delta/cp_partot_delta

    !write(*,*) 'IAM\CP:'
    !write(*,10) 'IAM\CP: HANS2:     ', cp_HANS2_delta/3600., ', ', cp_HANS2_delta/cp_HANS2_delta
    !write(*,10) 'IAM\CP: HANS2a:    ', cp_HANS2a_delta/3600., ', ', cp_HANS2a_delta/cp_HANS2_delta
    !write(*,10) 'IAM\CP: HANS2b:    ', cp_HANS2b_delta/3600., ', ', cp_HANS2b_delta/cp_HANS2_delta
    !write(*,10) 'IAM\CP: HANS2c:    ', cp_HANS2c_delta/3600., ', ', cp_HANS2c_delta/cp_HANS2_delta

    !write(*,*) 'IAM\CP:'
    !write(*,10) 'IAM\CP: mr2:     ', cp_mr2_delta/3600., ', ', cp_mr2_delta/cp_mr2_delta
    !write(*,10) 'IAM\CP: mr2a:    ', cp_mr2a_delta/3600., ', ', cp_mr2a_delta/cp_mr2_delta
    !write(*,10) 'IAM\CP: mr2b:    ', cp_mr2b_delta/3600., ', ', cp_mr2b_delta/cp_mr2_delta
    !write(*,10) 'IAM\CP: mr2c:    ', cp_mr2c_delta/3600., ', ', cp_mr2c_delta/cp_mr2_delta

    10 format(1a17, e24.16, 1a2, f5.2)

end subroutine
    
!
! ------------------
! --- kill agent ---
! ------------------
!
subroutine kill_agent(isp1,ia1)
    integer(4), intent(in) :: isp1,ia1
    a_sa(1,isp1,ia1) = 0
    naf_sp(isp1) = naf_sp(isp1) + 1
    iaf_sp(isp1,naf_sp(isp1)) = ia1
end subroutine

!
! ---------------------------
! --- get new agent index ---
! ---------------------------
!
subroutine get_new_ia(isp1,ia2,ierr)
    integer(4), intent(in) :: isp1
    integer(4), intent(out) :: ia2, ierr
    ierr = 0
    if (naf_sp(isp1).eq.0) then
        nat_sp(1,isp1)=nat_sp(1,isp1)+1
        if (nat_sp(1,isp1).gt.dimna_sp) then
            ierr = 1
            goto 10
        endif
        ia2=nat_sp(1,isp1)
    else
        ia2=iaf_sp(isp1,naf_sp(isp1))
        naf_sp(isp1)=naf_sp(isp1)-1
    endif
10 continue
end subroutine     

!
! ---------------------
! --- record change ---
! ---------------------
!
subroutine rec_c(isp1,ia1,int1,nt5,nt2,imr,icn1,iXSN,ierr,fpr,fpd)
    integer(4), intent(in) :: isp1,ia1,int1,nt5(5),nt2,imr
    real(4), intent(in) :: fpr,fpd
    integer(4), intent(out) :: iXSN,ierr
    integer(4), intent(inout) :: icn1
    integer(4) :: nt1, ig1, ip1, int2
    integer(4) :: nta(3), ntb(3), nta_o(3)
    integer(4) :: iaa1, iaa2, iaa1_o
    integer(4) :: cod1, cod2, cod1_o
    real(8) :: Loth
    real(8) :: fpi, rni, si, cfpi
    real(4) :: r
    !integer(4) :: pFlag
    !pFlag = 0
    ig1 = a_gid(1,isp1,ia1)
    nt1 = nt5(3)
    !if (isp1.eq.1) then
    !write(*,*) 'IAM\RC: isp1, ia1, int1 = ', isp1, ia1, int1
    !write(*,*) 'IAM\RC: nt1, nt2, imr = ', nt1, nt2, imr
    !write(*,*) 'IAM\RC: fpr, fpd = ', fpr, fpd
    !write(*,*) 'IAM\RC: ig1 = ', ig1
    !endif
    !
    ! check dim
    !
    ierr = 0
    if ((icn1.eq.-1).and.((a_cn(1,isp1,ia1)+1).gt.dimnc)) then
        ierr = 1
        goto 90
    endif
    !
    ! find aa
    !
    ip1 = g_ip(ig1,int1)
    !if (isp1.eq.1) then
    !write(*,*) 'IAM\RC: Protein. ip1 = ', ip1
    !endif
    if (ip1.eq.-9) then
        iXSN = 0
        goto 50
    endif
    !if (isp1.eq.1) then
    !write(*,*) 'IAM\RC: start, stop = ', g_pntstart(ig1,ip1), g_pntstop(ig1,ip1)
    !endif
    if (g_pother(ig1,ip1).eq.-1) then
        int2 = g_pntstart(ig1,ip1)+(int1-g_pntstart(ig1,ip1))/3*3.
        !if (isp1.eq.1) then
        !write(*,*) 'IAM\RC: + strand. int2 = ', int2
        !endif
        if (int1.eq.int2) then
            nta(1) = nt5(3)
            nta(2) = nt5(4)
            nta(3) = nt5(5)
            ntb(1) = nt2
            ntb(2) = nt5(4)
            ntb(3) = nt5(5)
            nta_o(1) = g_nt_o(ig1,int1)
            nta_o(2) = g_nt_o(ig1,int1+1)
            nta_o(3) = g_nt_o(ig1,int1+2)
        elseif (int1.eq.int2+1) then
            nta(1) = nt5(2)
            nta(2) = nt5(3)
            nta(3) = nt5(4)
            ntb(1) = nt5(2)
            ntb(2) = nt2
            ntb(3) = nt5(4)
            nta_o(1) = g_nt_o(ig1,int1-1)
            nta_o(2) = g_nt_o(ig1,int1)
            nta_o(3) = g_nt_o(ig1,int1+1)
        else
            nta(1) = nt5(1)
            nta(2) = nt5(2)
            nta(3) = nt5(3)
            ntb(1) = nt5(1)
            ntb(2) = nt5(2)
            ntb(3) = nt2
            nta_o(1) = g_nt_o(ig1,int1-2)
            nta_o(2) = g_nt_o(ig1,int1-1)
            nta_o(3) = g_nt_o(ig1,int1)
        endif
        cod1 = nta(1)*100+nta(2)*10+nta(3)*1
        cod2 = ntb(1)*100+ntb(2)*10+ntb(3)*1
        cod1_o = nta_o(1)*100+nta_o(2)*10+nta_o(3)*1
    else
        int2 = g_pntstop(ig1,ip1)-(g_pntstop(ig1,ip1)-int1)/3*3.
        !if (isp1.eq.1) then
        !write(*,*) 'IAM\RC: - strand. int2 = ', int2
        !endif
        if (int1.eq.int2) then
            nta(1) = nt5(3)
            nta(2) = nt5(2)
            nta(3) = nt5(1)
            ntb(1) = nt2
            ntb(2) = nt5(2)
            ntb(3) = nt5(1)
            nta_o(1) = g_nt_o(ig1,int1)
            nta_o(2) = g_nt_o(ig1,int1-1)
            nta_o(3) = g_nt_o(ig1,int1-2)
        elseif (int1.eq.int2-1) then
            nta(1) = nt5(4)
            nta(2) = nt5(3)
            nta(3) = nt5(2)
            ntb(1) = nt5(4)
            ntb(2) = nt2
            ntb(3) = nt5(2)
            nta_o(1) = g_nt_o(ig1,int1+1)
            nta_o(2) = g_nt_o(ig1,int1)
            nta_o(3) = g_nt_o(ig1,int1-1)
        else
            nta(1) = nt5(5)
            nta(2) = nt5(4)
            nta(3) = nt5(3)
            ntb(1) = nt5(5)
            ntb(2) = nt5(4)
            ntb(3) = nt2
            nta_o(1) = g_nt_o(ig1,int1+2)
            nta_o(2) = g_nt_o(ig1,int1+1)
            nta_o(3) = g_nt_o(ig1,int1)
        endif
        cod1 = int_opo(nta(1))*100+int_opo(nta(2))*10+int_opo(nta(3))*1
        cod2 = int_opo(ntb(1))*100+int_opo(ntb(2))*10+int_opo(ntb(3))*1
        cod1_o = int_opo(nta_o(1))*100+int_opo(nta_o(2))*10+int_opo(nta_o(3))*1
    endif ! pother
    !if (isp1.eq.1) then
    !write(*,*) 'IAM\RC: From codon = ', cod1
    !write(*,*) 'IAM\RC: To codon = ', cod2
    !write(*,*) 'IAM\RC: Original codon = ', cod1_o
    !endif
    iaa1 = gc(cod1)
    iaa2 = gc(cod2)
    iaa1_o = gc(cod1_o)
    !if (isp1.eq.1) then
    !write(*,*) 'IAM\RC: From aa = ', iaa1
    !write(*,*) 'IAM\RC: To aa = ', iaa2
    !write(*,*) 'IAM\RC: Original aa = ', iaa1_o
    !endif
    if (iaa1_o.eq.iaa2) then
        !if (isp1.eq.1) then
        !write(*,*) 'IAM\RC: Synonymous.'
        !endif
        iXSN = 1
    else
        !if (isp1.eq.1) then
        !write(*,*) 'IAM\RC: Nonsynonymous.'
        !endif
        iXSN = 2
    endif
    if (iXSN.eq.2) then
        if (oFTSS.eq.1) then
            if ((iaa1.gt.20).or.(iaa2.gt.20)) then
                !if (isp1.eq.1) then
                !write(*,*) 'IAM\RC: Change from/to START/STOP, killing agent...'
                !endif
                iXSN = -1
                if (oHANS.ge.0) then
                    call kill_agent(isp1,ia1)
                endif
                goto 90
            endif
        endif
    endif
    !if (isp1.eq.1) then
    !pause
    !endif

50 continue
    fpi = 1.
    cfpi = fpi
    if (oLM.gt.0) then
        if (iXSN.eq.0) then
            fpi = fpx
            cfpi = fpi
            rni = rnx
        elseif (iXSN.eq.1) then
            fpi = fps
            cfpi = fpi
            rni = rns
        elseif (iXSN.eq.2) then
            fpi = fpn
            cfpi = fpi
            rni = rnn
        else
            ierr = 2
            goto 80
        endif
        if (fpi.eq.1.) then
            goto 70
        endif
        if (oLM.eq.2) then
            r = r4_uni(idum_sp(isp1))
            if (r.eq.0.) then
                fpi = 0.
            else
                si = 1.-fpi
                si = -log(r)*si
                fpi = 1.-si
            endif
            cfpi = fpi
        elseif (oLM.eq.3) then
            !si = 1.-fpi
            !if (oKimura.eq.1) then
            !    si = random_gamma(fpbeta,rfirst) * si / fpbeta
            !else
            !    r = r4_uni(idum_sp(isp1))
            !    si = kimura(min(1+int(r*dimnk),dimnk))*si/fpbeta
            !endif
            !fpi = 1.-si
            !cfpi = fpi
        elseif (oLM.eq.4) then
            if (imr.eq.2) then ! assume N changes from reco are neutral
                fpi = 1.
                cfpi = fpi
                goto 70
            endif
            r = r4_uni(idum_sp(isp1))
            if (r.lt.rni) then
                fpi = 1.
                cfpi = fpi
                goto 70
            endif
        elseif (oLM.eq.5) then
            if (imr.eq.2) then ! assume N changes from reco are neutral
                fpi = 1.
                cfpi = fpi
                !if (isp1.eq.1) then
                !write(*,*) 'IAM\RC: Keeping due to reco.'
                !endif
                goto 70
            endif
            if ((iaa1.gt.20).or.(iaa2.gt.20)) then
                fpi = 1.
                cfpi = fpi
                !if (isp1.eq.1) then
                !write(*,*) 'IAM\RC: Keeping to STOP.'
                !endif
                goto 70
            endif
            if (MdXX(iaa1,iaa2).gt.Mdc) then
                rni = 0.
                !if (isp1.eq.1) then
                !write(*,*) 'IAM\RC: Killing MdXX > Mdc.'
                !endif
                goto 60
            endif
            rni = Mfc * (1. - MdXX(iaa1,iaa2) / Mdc)
            r = r4_uni(idum_sp(isp1))
            if (r.lt.rni) then
                fpi = 1.
                cfpi = fpi
                !if (isp1.eq.1) then
                !write(*,*) 'IAM\RC: Keeping by chance.'
                !endif
                goto 70
            endif
            !if (isp1.eq.1) then
            !write(*,*) 'IAM\RC: Killing by chance.'
            !endif
            elseif (oLM.eq.6) then
            if (imr.eq.2) then ! assume N changes from reco are neutral
                fpi = 1.
                cfpi = fpi
                !if (isp1.eq.1) then
                !write(*,*) 'IAM\RC: Keeping due to reco.'
                !endif
                goto 70
            endif
            if ((iaa1.gt.20).or.(iaa2.gt.20)) then
                fpi = 1.
                cfpi = fpi
                !if (isp1.eq.1) then
                !write(*,*) 'IAM\RC: Keeping to STOP.'
                !endif
                goto 70
            endif
            r = r4_uni(idum_sp(isp1))
            !if (isp1.eq.1) then
            !write(*,*) 'IAM\RC1:', r, fpi
            !endif
            si = 1.-fpi
            !if (isp1.eq.1) then
            !write(*,*) 'IAM\RC2:', si
            !endif
            si = si * MdXX(iaa1,iaa2) / Mdc
            !if (isp1.eq.1) then
            !write(*,*) 'IAM\RC3:', si, MdXX(iaa1,iaa2), Mdc
            !endif
            if (r.le.fdel) then
                if (r.eq.0.) then
                    fpi = 0.
                else
                    r = r / fdel
                    si = -log(r)*si
                    fpi = 1.-si
                endif
                cfpi = fpi
                !if (isp1.eq.1) then
                !write(*,*) 'IAM\RC4:', r, si, fpi
                !endif
            else
                si = (1.-fdel)*2./(1./si*fdel)
                !if (isp1.eq.1) then
                !write(*,*) 'IAM\RC5:', si
                !endif
                r = (r-fdel)/(1.-fdel)
                si = -(si-SQRT((1.-r)*si**2))
                fpi = 1.-si
                fpi = min(fpi, &
                g_Vmax0C(ig1) / a_VmaxC(1,isp1,ia1), &
                g_Vmax0N(ig1) / a_VmaxN(1,isp1,ia1), &
                g_Vmax0P(ig1) / a_VmaxP(1,isp1,ia1) )
                cfpi = fpi
                !if (isp1.eq.1) then
                !write(*,*) 'IAM\RC6:', r, si, fpi
                !pause
                !endif
            endif
            !if (isp1.eq.1) then
            !pause
            !endif
            elseif (oLM.eq.7) then
            if (imr.eq.2) then ! assume N changes from reco are neutral
                fpi = 1.
                cfpi = fpi
                !if (isp1.eq.1) then
                !write(*,*) 'IAM\RC: Keeping due to reco.'
                !endif
                goto 70
            endif
            if ((iaa1.gt.20).or.(iaa2.gt.20)) then
                fpi = 1.
                cfpi = fpi
                !if (isp1.eq.1) then
                !write(*,*) 'IAM\RC: Keeping to STOP.'
                !endif
                goto 70
            endif
            r = r4_uni(idum_sp(isp1))
            !if (isp1.eq.1) then
            !write(*,*) 'IAM\RC1:', r, fpi
            !endif
            si = 1.-fpi
            !if (isp1.eq.1) then
            !write(*,*) 'IAM\RC2:', si
            !endif
            si = si * MdXX(iaa1,iaa2) / Mdc
            !if (isp1.eq.1) then
            !write(*,*) 'IAM\RC3:', si, MdXX(iaa1,iaa2), Mdc
            !endif
            if (r.le.fdel) then
                if (r.eq.0.) then
                    fpi = 0.
                else
                    r = r / fdel
                    si = -log(r)*si
                    fpi = 1.-si
                endif
                cfpi = fpi
                !if (isp1.eq.1) then
                !write(*,*) 'IAM\RC4:', r, si, fpi
                !endif
            else
                si = (1.-fdel)*(1.-fpi)/fdel*2./3.
                !if (isp1.eq.1) then
                !write(*,*) 'IAM\RC5:', si
                !pause
                !endif
                r = (r-fdel)/(1.-fdel)
                si = -(r*si*2.)
                fpi = 1.-si
                fpi = min(fpi, &
                g_Vmax0C(ig1) / a_VmaxC(1,isp1,ia1), &
                g_Vmax0N(ig1) / a_VmaxN(1,isp1,ia1), &
                g_Vmax0P(ig1) / a_VmaxP(1,isp1,ia1) )
                cfpi = fpi
                !if (isp1.eq.1) then
                !write(*,*) 'IAM\RC6:', r, si, fpi
                !pause
                !endif
            endif
            !if (isp1.eq.1) then
            !pause
            !endif
            elseif (oLM.eq.8) then
            if (imr.eq.2) then ! for N changes from reco
                fpi = 1.
                !undo any previous fitness factor at this location in recipient
                fpi = fpi / fpr
                !adopt any previous fitness factor at this location from donor
                fpi = fpi * fpd
                !if ((fpr.ne.1.0).or.(fpd.ne.1.0)) then
                !if (fpd.ne.1.0) then
                !    pFlag = 1
                !endif
                !if (pFlag.eq.1) then
                !    write(*,*) 'IAM\RC1: Reco change to fp.'
                !    write(*,*) 'IAM\RC1: ia1, int1 = ', ia1, int1
                !    write(*,*) 'IAM\RC1: fpr, fpd, fpi = ', fpr, fpd, fpi
                !    write(*,*) 'IAM\RC1: V0/V = ', &
                !        g_Vmax0C(ig1) / a_VmaxC(1,isp1,ia1), &
                !        g_Vmax0N(ig1) / a_VmaxN(1,isp1,ia1), &
                !        g_Vmax0P(ig1) / a_VmaxP(1,isp1,ia1)
                !endif
                fpi = min(fpi, &
                g_Vmax0C(ig1) / a_VmaxC(1,isp1,ia1), &
                g_Vmax0N(ig1) / a_VmaxN(1,isp1,ia1), &
                g_Vmax0P(ig1) / a_VmaxP(1,isp1,ia1) )
                !if (pFlag.eq.1) then
                !    write(*,*) 'IAM\RC2: fpi = ', fpi
                !endif
                cfpi = fpd
                goto 70
            endif
            if ((iaa1.gt.20).or.(iaa2.gt.20)) then
                fpi = 1.
                cfpi = fpi
                !if (isp1.eq.1) then
                !write(*,*) 'IAM\RC: Keeping to STOP.'
                !endif
                goto 70
            endif
            r = r4_uni(idum_sp(isp1))
            !if (isp1.eq.1) then
            !write(*,*) 'IAM\RC1:', r, fpi
            !endif
            si = 1.-fpi
            !if (isp1.eq.1) then
            !write(*,*) 'IAM\RC2:', si
            !endif
            si = si * MdXX(iaa1,iaa2) / Mdc
            !if (isp1.eq.1) then
            !write(*,*) 'IAM\RC3:', si, MdXX(iaa1,iaa2), Mdc
            !endif
            if (r.le.fdel) then
                if (r.eq.0.) then
                    fpi = 0.
                else
                    r = r / fdel
                    si = -log(r)*si
                    fpi = 1.-si
                endif
                cfpi = fpi
                !if (isp1.eq.1) then
                !write(*,*) 'IAM\RC4:', r, si, fpi
                !endif
            else
                si = (1.-fdel)*(1.-fpi)/fdel*2./3.
                !if (isp1.eq.1) then
                !write(*,*) 'IAM\RC5:', si
                !pause
                !endif
                r = (r-fdel)/(1.-fdel)
                si = -(r*si*2.)
                fpi = 1.-si
                fpi = min(fpi, &
                g_Vmax0C(ig1) / a_VmaxC(1,isp1,ia1), &
                g_Vmax0N(ig1) / a_VmaxN(1,isp1,ia1), &
                g_Vmax0P(ig1) / a_VmaxP(1,isp1,ia1) )
                cfpi = fpi
                !if (isp1.eq.1) then
                !write(*,*) 'IAM\RC6:', r, si, fpi
                !pause
                !endif
            endif
            !if (isp1.eq.1) then
            !pause
            !endif
        endif
60 continue
        if (fpi.eq.0.) then
            iXSN = -2
            if (oHANS.ge.0) then
                call kill_agent(isp1,ia1)
            endif
            goto 90
        endif
70 continue    
    endif ! oLM>0
    !
    ! record change
    !
    if (icn1.eq.-1) then
        icn1 = a_cn(1,isp1,ia1) + 1
        a_cn(1,isp1,ia1) = icn1
        if (icn1.gt.ncg_sp(isp1)) then
            ncg_sp(isp1) = a_cn(1,isp1,ia1)
        endif
    endif
    a_ck(icn1,isp1,ia1,1) = int1
    a_cnt(icn1,isp1,ia1,1) = nt2
    a_cfp(icn1,isp1,ia1,1) = cfpi
    !if (pFlag.eq.1) then
    !    write(*,*) 'IAM\RC3: fpr, fpd, fpi = ', fpr, fpd, fpi
    !    write(*,*) 'IAM\RC3: a_cfp = ', a_cfp(icn1,isp1,ia1,1)
    !    write(*,*) 'IAM\RC3: V0/V = ', &
    !        g_Vmax0C(ig1) / a_VmaxC(1,isp1,ia1), &
    !        g_Vmax0N(ig1) / a_VmaxN(1,isp1,ia1), &
    !        g_Vmax0P(ig1) / a_VmaxP(1,isp1,ia1)
    !endif
        
    !
    ! update nnt & uEg
    !
    a_ntnk(1,isp1,ia1,nt1) = a_ntnk(1,isp1,ia1,nt1) - 1.
    a_ntnk(1,isp1,ia1,nt2) = a_ntnk(1,isp1,ia1,nt2) + 1.
    if (oSlim.eq.0) then
        if (g_ip(ig1,int1).ne.-9) then
            a_ntnSk(1,isp1,ia1,nt1) = a_ntnSk(1,isp1,ia1,nt1) - dble(g_ntfS(ig1,int1))/3.
            a_ntnSk(1,isp1,ia1,nt2) = a_ntnSk(1,isp1,ia1,nt2) + dble(g_ntfS(ig1,int1))/3.
            a_ntnNk(1,isp1,ia1,nt1) = a_ntnNk(1,isp1,ia1,nt1) - (1.-dble(g_ntfS(ig1,int1))/3.)
            a_ntnNk(1,isp1,ia1,nt2) = a_ntnNk(1,isp1,ia1,nt2) + (1.-dble(g_ntfS(ig1,int1))/3.)
        else
            a_ntnXk(1,isp1,ia1,nt1) = a_ntnXk(1,isp1,ia1,nt1) - 1.
            a_ntnXk(1,isp1,ia1,nt2) = a_ntnXk(1,isp1,ia1,nt2) + 1.
        endif
    endif
    a_uEg(1,isp1,ia1) = a_uEg(1,isp1,ia1) - uEbk(nt1) + uEbk(nt2)
    !
    ! update q0DNA
    !
    if ((oLQ.eq.11).or.(oLQ.eq.12).or.(oLQ.eq.31).or.(oLQ.eq.32).or.(oLQ.eq.33).or.(oLQ.eq.34)) then
        if (nt1.le.2) then
            a_q0DNAC(1,isp1,ia1) = a_q0DNAC(1,isp1,ia1) - rATC
            a_q0DNAN(1,isp1,ia1) = a_q0DNAN(1,isp1,ia1) - rATN
        else
            a_q0DNAC(1,isp1,ia1) = a_q0DNAC(1,isp1,ia1) - rGCC
            a_q0DNAN(1,isp1,ia1) = a_q0DNAN(1,isp1,ia1) - rGCN
        endif
        if (nt2.le.2) then
            a_q0DNAC(1,isp1,ia1) = a_q0DNAC(1,isp1,ia1) + rATC
            a_q0DNAN(1,isp1,ia1) = a_q0DNAN(1,isp1,ia1) + rATN
        else
            a_q0DNAC(1,isp1,ia1) = a_q0DNAC(1,isp1,ia1) + rGCC
            a_q0DNAN(1,isp1,ia1) = a_q0DNAN(1,isp1,ia1) + rGCN
        endif
    endif
    !
    ! update q0aa, nntt and Vmax
    !
    if (iXSN.eq.0) then
        goto 80
    endif
    if (iaa1.eq.iaa2) then
        goto 80
    endif
    if (iaa1.le.20) then
        a_naa(1,isp1,ia1,iaa1) = a_naa(1,isp1,ia1,iaa1) - 1
    endif
    if (iaa2.le.20) then
        a_naa(1,isp1,ia1,iaa2) = a_naa(1,isp1,ia1,iaa2) + 1
    endif
    if ((iaa1.le.20).and.(iaa2.le.20)) then
        a_naax(1,isp1,ia1,iaa1,iaa2) = a_naax(1,isp1,ia1,iaa1,iaa2) + 1
    endif
    if ((oLQ.eq.21).or.(oLQ.eq.22).or.(oLQ.eq.31).or.(oLQ.eq.32).or.(oLQ.eq.33).or.(oLQ.eq.34)) then
        if (iaa1.le.20) then
            a_q0aaC(1,isp1,ia1) = a_q0aaC(1,isp1,ia1) - m0aa / dble(g_naat(ig1)) * g_pgamma(ig1,ip1) * raaC(iaa1)
            a_q0aaN(1,isp1,ia1) = a_q0aaN(1,isp1,ia1) - m0aa / dble(g_naat(ig1)) * g_pgamma(ig1,ip1) * raaN(iaa1)
        endif
        if (iaa2.le.20) then
            a_q0aaC(1,isp1,ia1) = a_q0aaC(1,isp1,ia1) + m0aa / dble(g_naat(ig1)) * g_pgamma(ig1,ip1) * raaC(iaa2)
            a_q0aaN(1,isp1,ia1) = a_q0aaN(1,isp1,ia1) + m0aa / dble(g_naat(ig1)) * g_pgamma(ig1,ip1) * raaN(iaa2)
        endif
    endif
    if (oSlim.eq.0) then
        if (iaa1.le.20) then
            a_nntt(1,isp1,ia1,nta(1)) = a_nntt(1,isp1,ia1,nta(1)) - g_pgamma(ig1,ip1)
            a_nntt(1,isp1,ia1,nta(2)) = a_nntt(1,isp1,ia1,nta(2)) - g_pgamma(ig1,ip1)
            a_nntt(1,isp1,ia1,nta(3)) = a_nntt(1,isp1,ia1,nta(3)) - g_pgamma(ig1,ip1)
        endif
        if (iaa2.le.20) then
            a_nntt(1,isp1,ia1,ntb(1)) = a_nntt(1,isp1,ia1,ntb(1)) + g_pgamma(ig1,ip1)
            a_nntt(1,isp1,ia1,ntb(2)) = a_nntt(1,isp1,ia1,ntb(2)) + g_pgamma(ig1,ip1)
            a_nntt(1,isp1,ia1,ntb(3)) = a_nntt(1,isp1,ia1,ntb(3)) + g_pgamma(ig1,ip1)
        endif
    endif
    if (oLoth.eq.1) then
        Loth = 1.
        if (iaa1.le.20) then
            if (nta(1).ge.3) then
                Loth = Loth - g_pgamma(ig1,ip1) / g_nnttt(ig1) * sGCt
            endif
            if (nta(2).ge.3) then
                Loth = Loth - g_pgamma(ig1,ip1) / g_nnttt(ig1) * sGCt
            endif
            if (nta(3).ge.3) then
                Loth = Loth - g_pgamma(ig1,ip1) / g_nnttt(ig1) * sGCt
            endif
        endif
        if (iaa2.le.20) then
            if (ntb(1).ge.3) then
                Loth = Loth + g_pgamma(ig1,ip1) / g_nnttt(ig1) * sGCt
            endif
            if (ntb(2).ge.3) then
                Loth = Loth + g_pgamma(ig1,ip1) / g_nnttt(ig1) * sGCt
            endif
            if (ntb(3).ge.3) then
                Loth = Loth + g_pgamma(ig1,ip1) / g_nnttt(ig1) * sGCt
            endif
        endif
        a_VmaxC(1,isp1,ia1) = a_VmaxC(1,isp1,ia1) * Loth
        a_VmaxN(1,isp1,ia1) = a_VmaxN(1,isp1,ia1) * Loth
        a_VmaxP(1,isp1,ia1) = a_VmaxP(1,isp1,ia1) * Loth
    endif
80 continue
    if (oLM.gt.0) then
        a_VmaxC(1,isp1,ia1) = a_VmaxC(1,isp1,ia1) * fpi
        a_VmaxN(1,isp1,ia1) = a_VmaxN(1,isp1,ia1) * fpi
        a_VmaxP(1,isp1,ia1) = a_VmaxP(1,isp1,ia1) * fpi
        fpt_sp(isp1) = fpt_sp(isp1) + fpi * a_sr(1,isp1,ia1)
        fpn2_sp(isp1) = fpn2_sp(isp1) + a_sr(1,isp1,ia1)
    endif

    !if (pFlag.eq.1) then
    !    write(*,*) 'IAM\RC4: V0/V = ', &
    !            g_Vmax0C(ig1) / a_VmaxC(1,isp1,ia1), &
    !            g_Vmax0N(ig1) / a_VmaxN(1,isp1,ia1), &
    !            g_Vmax0P(ig1) / a_VmaxP(1,isp1,ia1)
    !    !if (fpd.ne.1.0) then
    !    !pause
    !    !endif
    !endif
    
    !
    ! XXX
    !
    if (g_nt_o(ig1,int1).ne.nt2) then
        a_cmr(icn1,isp1,ia1,1) = imr
        a_cXSN(icn1,isp1,ia1,1) = iXSN
    else ! un-mutated or un-recombined
        a_cmr(icn1,isp1,ia1,1) = 0
    endif

    
90 continue
end subroutine

!
! -------------------------
! --- calc uptake rates ---
! -------------------------
!
subroutine calc_V(isp1,ia1,LNC,LNN,LNP)
    integer(4), intent(in) :: isp1,ia1
    real(8), intent(in) :: LNC,LNN,LNP
    a_VC(1,isp1,ia1) = a_VmaxC(1,isp1,ia1) * LP * LNC
    a_VN(1,isp1,ia1) = a_VmaxN(1,isp1,ia1) * LP * LNN
    a_VP(1,isp1,ia1) = a_VmaxP(1,isp1,ia1) * LP * LNP
end subroutine

!
! -----------------------
! --- calc min quotas ---
! -----------------------
!
subroutine calc_q0(isp1,ia1,q0C,q0N,q0P)
    integer(4), intent(in) :: isp1,ia1
    real(8), intent(out) :: q0C, q0N, q0P
    q0C = q0othC
    q0N = q0othN
    q0P = q0othP
    ig1 = a_gid(1,isp1,ia1)
    if ((oLQ.eq.11).or.(oLQ.eq.31).or.(oLQ.eq.33)) then
        !q0C = q0C + a_q0DNAC(1,isp1,ia1)
        !q0N = q0N + a_q0DNAN(1,isp1,ia1)
        !q0P = q0P + a_q0DNAP(1,isp1,ia1)
        q0C = q0C + g_q0DNAC(ig1) + (a_q0DNAC(1,isp1,ia1) - g_q0DNAC(ig1)) * fdq0
        q0N = q0N + g_q0DNAN(ig1) + (a_q0DNAN(1,isp1,ia1) - g_q0DNAN(ig1)) * fdq0
        q0P = q0P + g_q0DNAP(ig1) + (a_q0DNAP(1,isp1,ia1) - g_q0DNAP(ig1)) * fdq0
    else
        q0C = q0C + g_q0DNAC(ig1)
        q0N = q0N + g_q0DNAN(ig1)
        q0P = q0P + g_q0DNAP(ig1)
    endif
    if ((oLQ.eq.21).or.(oLQ.eq.31).or.(oLQ.eq.34)) then
        !q0C = q0C + a_q0aaC(1,isp1,ia1)
        !q0N = q0N + a_q0aaN(1,isp1,ia1)
        q0C = q0C + g_q0aaC(ig1) + (a_q0aaC(1,isp1,ia1) - g_q0aaC(ig1)) * fdq0
        q0N = q0N + g_q0aaN(ig1) + (a_q0aaN(1,isp1,ia1) - g_q0aaN(ig1)) * fdq0
    else
        q0C = q0C + g_q0aaC(ig1)
        q0N = q0N + g_q0aaN(ig1)
    endif
end subroutine

!
! ------------------
! --- copy agent ---
! ------------------
!
subroutine copy_agent(ix1,isp1,ia1,ix2,isp2,ia2)
    integer(4), intent(in) :: ix1,isp1,ia1,ix2,isp2,ia2
    a_sa(ix2,isp2,ia2) = a_sa(ix1,isp1,ia1)
    a_sr(ix2,isp2,ia2) = a_sr(ix1,isp1,ia1)
    a_srm(ix2,isp2,ia2) = a_srm(ix1,isp1,ia1)
    a_srr(ix2,isp2,ia2) = a_srr(ix1,isp1,ia1)
    a_iar(ix2,isp2,ia2) = a_iar(ix1,isp1,ia1)
    a_gid(ix2,isp2,ia2) = a_gid(ix1,isp1,ia1)
    a_cn(ix2,isp2,ia2) = a_cn(ix1,isp1,ia1)
    do ic1 = 1, a_cn(ix1,isp1,ia1)
        a_ck(ic1,isp2,ia2,ix2) = a_ck(ic1,isp1,ia1,ix1)
    enddo
    do ic1 = 1, a_cn(ix1,isp1,ia1)
        a_cnt(ic1,isp2,ia2,ix2) = a_cnt(ic1,isp1,ia1,ix1)
        a_cfp(ic1,isp2,ia2,ix2) = a_cfp(ic1,isp1,ia1,ix1)
    enddo
    do nt1 = 1, 4
        a_ntnk(ix2,isp2,ia2,nt1) = a_ntnk(ix1,isp1,ia1,nt1)
    enddo
    do iaa1 = 1, 20
        a_naa(ix2,isp2,ia2,iaa1) = a_naa(ix1,isp1,ia1,iaa1)
        do iaa2 = 1, 20
            a_naax(ix2,isp2,ia2,iaa1,iaa2) = a_naax(ix1,isp1,ia1,iaa1,iaa2)
        enddo
    enddo
    a_uEg(ix2,isp2,ia2) = a_uEg(ix1,isp1,ia1)
    a_q0DNAC(ix2,isp2,ia2) = a_q0DNAC(ix1,isp1,ia1)
    a_q0DNAN(ix2,isp2,ia2) = a_q0DNAN(ix1,isp1,ia1)
    a_q0DNAP(ix2,isp2,ia2) = a_q0DNAP(ix1,isp1,ia1)
    a_q0aaC(ix2,isp2,ia2) = a_q0aaC(ix1,isp1,ia1)
    a_q0aaN(ix2,isp2,ia2) = a_q0aaN(ix1,isp1,ia1)
    a_VmaxC(ix2,isp2,ia2) = a_VmaxC(ix1,isp1,ia1)
    a_VmaxN(ix2,isp2,ia2) = a_VmaxN(ix1,isp1,ia1)
    a_VmaxP(ix2,isp2,ia2) = a_VmaxP(ix1,isp1,ia1)
    if (oSlim.eq.0) then
        a_sn(ix2,isp2,ia2) = a_sn(ix1,isp1,ia1)
        a_rsn(ix2,isp2,ia2) = a_rsn(ix1,isp1,ia1)
        a_nd(ix2,isp2,ia2) = a_nd(ix1,isp1,ia1)
        do nt1 = 1, 4
            a_ntnSk(ix2,isp2,ia2,nt1) = a_ntnSk(ix1,isp1,ia1,nt1)
            a_ntnNk(ix2,isp2,ia2,nt1) = a_ntnNk(ix1,isp1,ia1,nt1)
            a_ntnXk(ix2,isp2,ia2,nt1) = a_ntnXk(ix1,isp1,ia1,nt1)
        enddo
        a_nmt(ix2,isp2,ia2) = a_nmt(ix1,isp1,ia1)
        a_nmS(ix2,isp2,ia2) = a_nmS(ix1,isp1,ia1)
        a_nmN(ix2,isp2,ia2) = a_nmN(ix1,isp1,ia1)
        a_nrt(ix2,isp2,ia2) = a_nrt(ix1,isp1,ia1)
        a_nrS(ix2,isp2,ia2) = a_nrS(ix1,isp1,ia1)
        a_nrN(ix2,isp2,ia2) = a_nrN(ix1,isp1,ia1)
        do ic1 = 1, a_cn(ix1,isp1,ia1)
            a_cXSN(ic1,isp2,ia2,ix2) = a_cXSN(ic1,isp1,ia1,ix1)
        enddo
        do ic1 = 1, a_cn(ix1,isp1,ia1)
            a_cmr(ic1,isp2,ia2,ix2) = a_cmr(ic1,isp1,ia1,ix1)
        enddo
        do nt1 = 1, 4
            do nt2 = 1, 4
                a_nm(ix2,isp2,ia2,nt1,nt2) = a_nm(ix1,isp1,ia1,nt1,nt2)
                a_nr(ix2,isp2,ia2,nt1,nt2) = a_nr(ix1,isp1,ia1,nt1,nt2)
            enddo
        enddo
        do nt1 = 1, 4
            a_nntt(ix2,isp2,ia2,nt1) = a_nntt(ix1,isp1,ia1,nt1)
        enddo
        a_VC(ix2,isp2,ia2) = a_VC(ix1,isp1,ia1)
        a_VN(ix2,isp2,ia2) = a_VN(ix1,isp1,ia1)
        a_VP(ix2,isp2,ia2) = a_VP(ix1,isp1,ia1)
        a_qC(ix2,isp2,ia2) = a_qC(ix1,isp1,ia1)
        a_qN(ix2,isp2,ia2) = a_qN(ix1,isp1,ia1)
        a_qP(ix2,isp2,ia2) = a_qP(ix1,isp1,ia1)
        a_CLimit(ix2,isp2,ia2) = a_CLimit(ix1,isp1,ia1)
        a_NLimit(ix2,isp2,ia2) = a_NLimit(ix1,isp1,ia1)
        a_PLimit(ix2,isp2,ia2) = a_PLimit(ix1,isp1,ia1)
        a_tb(ix2,isp2,ia2) = a_tb(ix1,isp1,ia1)
        a_tg(ix2,isp2,ia2) = a_tg(ix1,isp1,ia1)
    endif
end subroutine     

!
! -----------------------------
! --- calc kg & other stuff ---
! -----------------------------
!
subroutine calc_A(isp1,ia1,q0C,q0N,q0P,kg)
    integer(4), intent(in) :: isp1,ia1
    real(8), intent(in) :: q0C,q0N,q0P
    real(8), intent(out) :: kg
    qaveC = q0C * log(2.) * 2.
    qaveN = q0N * log(2.) * 2.
    qaveP = q0P * log(2.) * 2.
    a_qC(1,isp1,ia1) = qaveC
    a_qN(1,isp1,ia1) = qaveN
    a_qP(1,isp1,ia1) = qaveP
    kgC = a_VC(1,isp1,ia1) / qaveC
    kgN = a_VN(1,isp1,ia1) / qaveN
    kgP = a_VP(1,isp1,ia1) / qaveP
    a_CLimit(1,isp1,ia1) = 2.
    a_NLimit(1,isp1,ia1) = 2.
    a_PLimit(1,isp1,ia1) = 2.
    if ((kgC.le.kgN).and.(kgC.le.kgP)) then
        kg = kgC
        a_CLimit(1,isp1,ia1) = 1.
    elseif ((kgN.le.kgC).and.(kgN.le.kgP)) then
        kg = kgN
        a_NLimit(1,isp1,ia1) = 1.
    elseif ((kgP.le.kgC).and.(kgP.le.kgN)) then
        kg = kgP
        a_PLimit(1,isp1,ia1) = 1.
    else
        write(*,*) 'IAM: ERR in calc_A.'
        write(*,*) 'IAM: isp1, ia1 = ', isp1, ia1
        write(*,*) 'IAM: kg = ', kgC, kgN, kgP
        write(*,*) 'IAM: V = ', a_VC(1,isp1,ia1), a_VN(1,isp1,ia1), a_VP(1,isp1,ia1)
        write(*,*) 'IAM: qave = ', qaveC, qaveN, qaveP
    endif
end subroutine

subroutine calc_A_Slim(isp1,ia1,LNC,LNN,LNP,kg)
    integer(4), intent(in) :: isp1,ia1
    real(8), intent(in) :: LNC,LNN,LNP
    real(8), intent(out) :: kg
    real(4) VC, VN, VP
    real(8) q0C, q0N, q0P
    VC = a_VmaxC(1,isp1,ia1) * LP * LNC
    VN = a_VmaxN(1,isp1,ia1) * LP * LNN
    VP = a_VmaxP(1,isp1,ia1) * LP * LNP
    q0C = q0othC
    q0N = q0othN
    q0P = q0othP
    ig1 = a_gid(1,isp1,ia1)
    if ((oLQ.eq.11).or.(oLQ.eq.31).or.(oLQ.eq.33)) then
        !q0C = q0C + a_q0DNAC(1,isp1,ia1)
        !q0N = q0N + a_q0DNAN(1,isp1,ia1)
        !q0P = q0P + a_q0DNAP(1,isp1,ia1)
        q0C = q0C + g_q0DNAC(ig1) + (a_q0DNAC(1,isp1,ia1) - g_q0DNAC(ig1)) * fdq0
        q0N = q0N + g_q0DNAN(ig1) + (a_q0DNAN(1,isp1,ia1) - g_q0DNAN(ig1)) * fdq0
        q0P = q0P + g_q0DNAP(ig1) + (a_q0DNAP(1,isp1,ia1) - g_q0DNAP(ig1)) * fdq0
    else
        q0C = q0C + g_q0DNAC(ig1)
        q0N = q0N + g_q0DNAN(ig1)
        q0P = q0P + g_q0DNAP(ig1)
    endif
    if ((oLQ.eq.21).or.(oLQ.eq.31).or.(oLQ.eq.34)) then
        !q0C = q0C + a_q0aaC(1,isp1,ia1)
        !q0N = q0N + a_q0aaN(1,isp1,ia1)
        q0C = q0C + g_q0aaC(ig1) + (a_q0aaC(1,isp1,ia1) - g_q0aaC(ig1)) * fdq0
        q0N = q0N + g_q0aaN(ig1) + (a_q0aaN(1,isp1,ia1) - g_q0aaN(ig1)) * fdq0
    else
        q0C = q0C + g_q0aaC(ig1)
        q0N = q0N + g_q0aaN(ig1)
    endif
    qaveC = q0C * log(2.) * 2.
    qaveN = q0N * log(2.) * 2.
    qaveP = q0P * log(2.) * 2.
    kgC = VC / qaveC
    kgN = VN / qaveN
    kgP = VP / qaveP
    if ((kgC.le.kgN).and.(kgC.le.kgP)) then
        kg = kgC
    elseif ((kgN.le.kgC).and.(kgN.le.kgP)) then
        kg = kgN
    elseif ((kgP.le.kgC).and.(kgP.le.kgN)) then
        kg = kgP
    else
        write(*,*) 'IAM: ERR in calc_A_Slim.'
        write(*,*) 'IAM: isp1, ia1 = ', isp1, ia1
        write(*,*) 'IAM: kg = ', kgC, kgN, kgP
        write(*,*) 'IAM: V = ', VC, VN, VP
        write(*,*) 'IAM: qave = ', qaveC, qaveN, qaveP
    endif
end subroutine

!
! ------------------
! --- get new sn ---
! ------------------
!
subroutine new_sn(ix1,isp1,ia1)
    integer(4), intent(in) :: ix1,isp1,ia1
    if (oSlim.eq.0) then
        a_sn(ix1,isp1,ia1) = gsn_sp(isp1)
        gsn_sp(isp1) = gsn_sp(isp1) + 1
        if (gsn_sp(isp1).gt.gsnmax) then
            !write(*,*) 'IAM\NEWSN: Warning: resetting sn. isp1 = ', isp1
            gsn_sp(isp1) = 1
        endif
    endif
end subroutine     

!
! ---------------------------
! --- fixed changes count ---
! ---------------------------
!
subroutine fix_changes()
    cp_fix_t1 = omp_get_wtime()
    do ig1 = 1, ngt
        ifx(ig1) = 0
        do isp1 = 1, dimnsp
            ifx_sp(isp1,ig1) = 0
        enddo
    enddo
    do isp1 = 1, dimnsp
        do ia1 = 1, nat_sp(1,isp1)
            if (a_sa(1,isp1,ia1).ge.1) then
                ig1 = a_gid(1,isp1,ia1)
                !global
                if (ifx(ig1).eq.0) then
                    ifx(ig1) = 1
                    do ic1 = 1, a_cn(1,isp1,ia1)
                        fx_ck(ic1,ig1) = a_ck(ic1,isp1,ia1,1)
                        fx_cnt(ic1,ig1) = a_cnt(ic1,isp1,ia1,1)
                        if (oSlim.eq.0) then
                            fx_cXSN(ic1,ig1) = a_cXSN(ic1,isp1,ia1,1)
                            fx_cmr(ic1,ig1) = a_cmr(ic1,isp1,ia1,1)
                        endif
                    enddo
                    fx_cn(ig1) = a_cn(1,isp1,ia1)
                else
                    fx_cn(ig1) = min(fx_cn(ig1), a_cn(1,isp1,ia1))
                    do ic1 = 1, fx_cn(ig1)
                        if (a_ck(ic1,isp1,ia1,1).ne.fx_ck(ic1,ig1)) then
                        !if ((a_ck(ic1,isp1,ia1,1).ne.fx_ck(ic1,ig1)).or.(a_cnt(ic1,isp1,ia1,1).ne.fx_cnt(ic1,ig1))) then
                            fx_cn(ig1) = ic1 - 1
                            goto 10
                        endif
                    enddo
                endif
10 continue
                !by sp
                if (ifx_sp(isp1,ig1).eq.0) then
                    ifx_sp(isp1,ig1) = 1
                    do ic1 = 1, a_cn(1,isp1,ia1)
                        fx_ck_sp(ic1,isp1,ig1) = a_ck(ic1,isp1,ia1,1)
                        fx_cnt_sp(ic1,isp1,ig1) = a_cnt(ic1,isp1,ia1,1)
                        if (oSlim.eq.0) then
                            fx_cXSN_sp(ic1,isp1,ig1) = a_cXSN(ic1,isp1,ia1,1)
                            fx_cmr_sp(ic1,isp1,ig1) = a_cmr(ic1,isp1,ia1,1)
                        endif
                    enddo
                    fx_cn_sp(isp1,ig1) = a_cn(1,isp1,ia1)
                else
                    fx_cn_sp(isp1,ig1) = min(fx_cn_sp(isp1,ig1), a_cn(1,isp1,ia1))
                    do ic1 = 1, fx_cn_sp(isp1,ig1)
                        if (a_ck(ic1,isp1,ia1,1).ne.fx_ck_sp(ic1,isp1,ig1)) then
                        !if ((a_ck(ic1,isp1,ia1,1).ne.fx_ck_sp(ic1,isp1,ig1)).or.(a_cnt(ic1,isp1,ia1,1).ne.fx_cnt_sp(ic1,isp1,ig1))) then
                            fx_cn_sp(isp1,ig1) = ic1 - 1
                            goto 20
                        endif
                    enddo
                endif
20 continue
            endif ! sa
        enddo ! ia1
    enddo ! isp1
    g_ncft = 0
    do ig1 = 1, ngt
        g_ncf(ig1) = 0
        g_ncSf(ig1) = 0
        g_ncNf(ig1) = 0
        g_ncmf(ig1) = 0
        g_ncrf(ig1) = 0
        do ic1 = 1, fx_cn(ig1)
            g_ncf(ig1) = g_ncf(ig1) + 1
            g_ncft = g_ncft + 1
            if (fx_cXSN(ic1,ig1).eq.1) then
                g_ncSf(ig1) = g_ncSf(ig1) + 1
            endif
            if (fx_cXSN(ic1,ig1).eq.2) then
                g_ncNf(ig1) = g_ncNf(ig1) + 1
            endif
            if (fx_cmr(ic1,ig1).eq.1) then
                g_ncmf(ig1) = g_ncmf(ig1) + 1
            endif
            if (fx_cmr(ic1,ig1).eq.2) then
                g_ncrf(ig1) = g_ncrf(ig1) + 1
            endif
        enddo
        write(*,*) 'IAM\FIX: Changes fixed. ig1 = ', ig1
        write(*,*) 'IAM\FIX: ncf/f  = ', g_ncf(ig1), g_ncff(ig1)
        write(*,*) 'IAM\FIX: ncSf/f = ', g_ncSf(ig1), g_ncSff(ig1)
        write(*,*) 'IAM\FIX: ncNf/f = ', g_ncNf(ig1), g_ncNff(ig1)
        write(*,*) 'IAM\FIX: ncmf/f = ', g_ncmf(ig1), g_ncmff(ig1)
        write(*,*) 'IAM\FIX: ncrf/f = ', g_ncrf(ig1), g_ncrff(ig1)
        write(*,*) 'IAM\FIX: isp1, fx_cn_sp:'
        do isp1 = 1, dimnsp
            write(*,*) 'IAM\FIX: ', isp1, fx_cn_sp(isp1,ig1)
        enddo
    enddo ! ig1
    cp_fix_delta = cp_fix_delta + omp_get_wtime() - cp_fix_t1
end subroutine

end
    
!
! -----------------
! -----------------
! --- functions ---
! -----------------
! -----------------
!
! --------------
! --- sarkar ---
! --------------
!
function sarkar(mx, nx, idum)
    real(8) mx, nx
    real(4) r, CDFx
    real(8), allocatable :: px(:)
    integer(4) idum, sarkar, rx, ix, inx
    !real(8) erx
    !write(*,*) 'IAM\SARKAR\IN: mx, nx = ', mx, nx
    inx = int(nx)
    if (nx.gt.(2147483647.-1.)) then
        write(*,*) 'IAM\SARKAR: ERR: nx exceeds int(4) limit. nx = ', nx
        pause
    endif
    if (inx.eq.0) then
        rx = 0
        goto 10
    endif
    allocate(px(inx)) ! note: this takes a lot of time, should pre-allocate, but oHANSs = 1 takes care of this now
    r = r4_uni(idum)
    !write(*,*) 'IAM\SARKAR: r = ', r
    px(0+1) = exp(-mx)
    CDFx = px(0+1)
    !write(*,*) 'IAM\SARKAR: rx, px, CDFx = ', 0, px(0+1), CDFx
    !erx = 0. !dble(0)*px(0+1)
    if (r.le.CDFx) then
        rx = 0
        goto 10
    endif
    do rx = 1, inx
        if (rx.ge.inx) then
            goto 10
        endif
        px(rx+1) = 0.
        do ix = 0, rx - 1
            px(rx+1) = px(rx+1) + px(ix+1) / dble(rx - ix + 1)
        enddo
        px(rx+1) = mx / dble(rx) * px(rx+1)
        !erx = erx + dble(rx)*px(rx+1)
        CDFx = CDFx + px(rx+1)
        !write(*,*) 'IAM\SARKAR: rx, px, CDFx = ', rx, px(rx+1), CDFx
        if (r.le.CDFx) then
            !pause
            goto 10
        endif
    enddo
    write(*,*) 'IAM\SARKAR: ERR: mx, nx, CDFx, r ', mx, nx, CDFx, r, rx, inx
    pause
10 continue
sarkar = rx
!sarkar = erx*1.e6
!write(*,*) 'IAM\SARKAR\OUT: rx = ', rx
!write(*,*) 'IAM\SARKAR\OUT: erx, CDFx = ', erx, CDFx
!pause
return
end

!
! -------------
! --- bucci ---
! -------------
!
function bucci(hx, nx, idum)
    real(8) bucci, hx, nx, px
    real(8) CDFx
    real(4) r, r1, r2, n1, n2
    integer(4) idum
    integer(4) ix
    if ((hx.le.0.).or.(nx.le.0.)) then
        x = 0.
        goto 10
    endif
    if (hx.le.(200.)) then
        r = r4_uni(idum)
        px = exp(-hx)
        CDFx = px
        x = 0.
        do ix = 1, nx
            if (r.le.CDFx) then
                goto 10
            endif
            x = x + 1.
            px = px * hx / dble(x)
            CDFx = CDFx + px
        enddo
    else
5 continue
        r1 = r4_uni(idum)
        if (r1.eq.0.) goto 5
6 continue
        r2 = r4_uni(idum)
        if (r2.eq.0.) goto 6
        n1 = sqrt(-2. * log(r1)) * cos(2. * 3.14159265 * r2)
        !n2 = sqrt(-2. * log(r1)) * sin(2. * 3.14159265 * r2)
        x = hx + (n1 * sqrt(hx))
    endif
10 continue
bucci = x
return
end

!
! -----------    
! --- sap ---
! -----------    
!    
function sap(iap)
    character(4) sap
    character(1) sap1
    character(2) sap2
    character(3) sap3
    character(4) sap4
    integer(4) iap
    if (iap.le.9) then
        write(sap1,'(i1)') iap
        sap4='000'//sap1
    elseif (iap.le.99) then
        write(sap2,'(i2)') iap
        sap4='00'//sap2
    elseif (iap.le.999) then
        write(sap3,'(i3)') iap
        sap4='0'//sap3
    elseif (iap.le.9999) then
        write(sap4,'(i4)') iap
    else
        write(*,*) 'IAM\SAP: ERR: iap > 9999. iap = ', iap
        sap4 = "XXXX"
    endif
    sap = sap4
end function

!
! -------------
! --- kreft ---
! -------------
!
function kreft_p(Mean, CV, CVLimit, LLimit, HLimit, idum,kn,fn,wn)
real(8) Mean, CV, CVLimit, LLimit, HLimit
real(8) kreft_p
integer(4) k, kMax
integer(4), INTENT(INOUT) :: idum
real(4), INTENT(IN) :: fn(128)
integer(4), INTENT(IN) :: kn(128)
real(4), INTENT(IN) :: wn(128)
if (CV.eq.0.0) then
    goto 20
endif
kMax = 100
do 10 k = 1, kMax
    kreft_p = r4_nor(idum,kn,fn,wn) * CV * Mean + Mean
    if ((kreft_p.ge.(Mean - CVLimit * CV * Mean)) .and. &
        (kreft_p.le.(Mean + CVLimit * CV * Mean)) .and. &
        (kreft_p.gt.LLimit) .and. &
        (kreft_p.le.HLimit)) then
    goto 30
    endif
10 continue
write(*,*) 'IAM: kreft: ERR: Problem finding valid value. Assigning mean...'
write(*,*) 'Mean, CV, CVLimit, LLimit, HLimit = '
write(*,*) Mean, CV, CVLimit, LLimit, HLimit
20 continue
kreft_p = Mean
30 continue
return
end

!
! --------------
! --- pcheck ---
! --------------
!
function pcheck(px,jdum)
integer(4) pcheck, jdum
real(8) px
real(8) ranprec
ranprec = 1e-3
pcheck = 0
if (px.gt.ranprec) then
    r = r4_uni(jdum)
    if (r.le.px) then
        pcheck = 1
    endif
else
    r = r4_uni(jdum)
    if (r.lt.ranprec) then
        r = r4_uni(jdum)
        if (r*ranprec.lt.px) then
            pcheck = 1
        endif
    endif
endif
end function

!
! XXXXXXXXXXX
! XXXXXXXXXXX
! XXX EOF XXX
! XXXXXXXXXXX
! XXXXXXXXXXX
!

