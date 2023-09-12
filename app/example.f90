module app_example
  use iso_fortran_env,only:wp => real64

  real(wp),parameter,private :: autoaa = 0.52917721067d0
  real(wp),parameter,private :: aatoau = 1.0d0/autoaa

!> prophyrine+AlCl
!&<
  integer,parameter :: testnat_porph = 58
  integer,parameter :: testat_porph(testnat_porph) = &
  & [17,13,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,&
  &   6, 6,6,6,6,6,6,6,6,6,6,7,7,7,7,7,7,7,7,1,1,1,1,&
  &   1, 1,1,1,1,1,1,1,1,1,1,1 ]
  real(wp),parameter :: testxyz_porph(3,testnat_porph) = reshape( &
  &    [17.80488069531964_wp, 16.97314082360259_wp, 3.95395834052218_wp,&
  &     17.77480726394300_wp, 16.99670975911163_wp, 8.01004862153323_wp,&
  &     17.55651926786093_wp,  4.56976712104640_wp, 9.17090143703508_wp,&
  &     20.21206206482194_wp,  4.81168263416407_wp, 9.19808029849216_wp,&
  &     15.99361266718554_wp,  6.68984114693105_wp, 9.05264765256453_wp,&
  &     21.37340819346931_wp,  7.17791906928089_wp, 9.10065952216448_wp,&
  &     17.15670059469320_wp,  9.06182563922934_wp, 8.95464616323316_wp,&
  &     19.80718143753246_wp,  9.30345011554356_wp, 8.97148028298884_wp,&
  &     16.14042277267272_wp, 11.61191133266435_wp, 8.87034515441769_wp,&
  &     20.34863926263168_wp, 11.99578050492211_wp, 8.88967482830795_wp,&
  &      7.96042046871797_wp, 13.38930160729190_wp, 9.18654143851209_wp,&
  &     12.77042742555772_wp, 14.42017658513013_wp, 8.89354925929706_wp,&
  &      5.59329930137183_wp, 14.54473453554468_wp, 9.30927799477929_wp,&
  &     10.08173990006294_wp, 14.95923836910088_wp, 9.02417551887075_wp,&
  &      9.83737562244289_wp, 17.60942286500477_wp, 9.00346655238609_wp,&
  &      5.34839924844815_wp, 17.20052288358834_wp, 9.27914335515604_wp,&
  &     12.38440772545960_wp, 18.62841595967898_wp, 8.87012413198713_wp,&
  &     22.76742546247186_wp, 19.57880177194378_wp, 8.93445035202396_wp,&
  &      7.46399488670201_wp, 18.76673930272034_wp, 9.13118820759162_wp,&
  &     15.18934324174975_wp, 22.00235509075167_wp, 8.84837509774140_wp,&
  &     19.39729791007948_wp, 22.38615114948483_wp, 8.88165877787674_wp,&
  &     15.72952892175636_wp, 24.69397036976622_wp, 8.94709761228124_wp,&
  &     18.37983141219312_wp, 24.93550570767372_wp, 8.96911211532866_wp,&
  &     14.16116923133441_wp, 26.81906725132794_wp, 9.05939044296441_wp,&
  &     19.54110677057694_wp, 27.30688811415975_wp, 9.09814760299700_wp,&
  &     15.32070241291261_wp, 29.18484351910493_wp, 9.18200320323787_wp,&
  &     17.97639886354862_wp, 29.42631286775734_wp, 9.19990114461879_wp,&
  &     23.15375471869472_wp, 15.37042352138243_wp, 8.93365498496856_wp,&
  &     28.07167918427908_wp, 15.23563744429889_wp, 9.25139840717964_wp,&
  &     25.69974255875718_wp, 16.39073589093223_wp, 9.08502834571098_wp,&
  &     25.45490271252499_wp, 19.04094480210461_wp, 9.08777179418927_wp,&
  &     30.18422134549823_wp, 16.80406633651126_wp, 9.41859376897403_wp,&
  &     27.57304736528073_wp, 20.61319687591900_wp, 9.26508333524511_wp,&
  &     29.93878759597024_wp, 19.45991455109562_wp, 9.42645384136645_wp,&
  &     13.67962246889538_wp, 12.08427678251964_wp, 8.88076008056737_wp,&
  &     22.68304640794519_wp, 12.90786843632262_wp, 8.91383505928998_wp,&
  &     18.10371190315771_wp, 13.34771382569003_wp, 8.83781948144163_wp,&
  &     14.12050191669306_wp, 16.66553342278291_wp, 8.82140575027951_wp,&
  &     21.41838941393265_wp, 17.33329643956978_wp, 8.86132681186142_wp,&
  &     17.43476992370002_wp, 20.64999954835533_wp, 8.81789297517925_wp,&
  &     12.85479289294967_wp, 21.09096057565307_wp, 8.84713413313874_wp,&
  &     21.85800692116170_wp, 21.91462239500707_wp, 8.91300577164098_wp,&
  &     13.94007937532288_wp,  6.52265389206403_wp, 9.03105635380341_wp,&
  &     23.42333184284292_wp,  7.38147400040357_wp, 9.11942284497367_wp,&
  &      8.17416635558616_wp, 11.34015930320184_wp, 9.21029033269129_wp,&
  &      3.89568924246263_wp, 13.37922516898219_wp, 9.42463065197876_wp,&
  &      3.46610633499495_wp, 18.03845754387010_wp, 9.36580484044416_wp,&
  &      7.29357629372672_wp, 20.81970424039168_wp, 9.10510746629872_wp,&
  &     12.11119344440215_wp, 26.61565518354103_wp, 9.04673201442255_wp,&
  &     21.59460083171406_wp, 27.47439853118350_wp, 9.11064527249161_wp,&
  &     14.15604544535623_wp, 30.88430966754304_wp, 9.26054827471168_wp,&
  &     18.81422061266223_wp, 31.30846688630078_wp, 9.28770149257794_wp,&
  &     32.06510905820445_wp, 15.96742336109581_wp, 9.54031174367235_wp,&
  &     28.24334768026415_wp, 13.18265498886139_wp, 9.24341100500465_wp,&
  &     27.35735196079739_wp, 22.66222581530651_wp, 9.27420542379521_wp,&
  &     31.63418593768974_wp, 20.62681212976239_wp, 9.55886268010183_wp,&
  &     16.71745260461593_wp,  2.68724590971498_wp, 9.23496326541442_wp,&
  &     21.37541873035842_wp,  3.11200077342733_wp, 9.28927837499231_wp],&
  &shape(testxyz_porph))
!&>

!> definitions
  integer,parameter :: testnat = testnat_porph
  integer,parameter :: testat(testnat) = testat_porph
  real(wp),parameter :: testxyz(3,testnat) = testxyz_porph

!&<
  character(len=2),dimension(118),parameter :: testelem = [&
  &    'H ','He','Li','Be','B ','C ','N ','O ','F ','Ne','Na',&
  &    'Mg','Al','Si','P ','S ','Cl','Ar','K ','Ca','Sc','Ti','V ','Cr','Mn',&
  &    'Fe','Co','Ni','Cu','Zn','Ga','Ge','As','Se','Br','Kr','Rb','Sr',&
  &    'Y ','Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In','Sn','Sb',&
  &    'Te','I ','Xe','Cs','Ba','La','Ce','Pr','Nd','Pm','Sm','Eu','Gd',&
  &    'Tb','Dy','Ho','Er','Tm','Yb','Lu','Hf','Ta','W ','Re','Os','Ir',&
  &    'Pt','Au','Hg','Tl','Pb','Bi','Po','At','Rn','Fr','Ra','Ac','Th',&
  &    'Pa','U ','Np','Pu','Am','Cm','Bk','Cf','Es','Fm','Md','No','Lr',&
  &    'Rf','Db','Sg','Bh','Hs','Mt','Ds','Rg','Cn','Nh','Fl','Mc','Lv',&
  &    'Ts','Og']
  private :: testelem
!&>

  public :: lwoniom_write_example

!========================================================================================!
!========================================================================================!
contains !> MODULE PROCEDURES START HERE
!========================================================================================!
!========================================================================================!

  subroutine lwoniom_write_example
    use iso_fortran_env,only:stdout => output_unit
    call write_example_xyz
    write (stdout,'(a,i0,a)') 'A test system with ',testnat,' atoms was written to example.xyz'
    call write_example_toml
    write (stdout,'(a)') 'A matching TOML input file was written to example.toml'
    write (stdout,'(a)') 'To run the example use:'
    write (stdout,'(a)') '  lwoniom-app example.toml'
    write (stdout,*)
    write (stdout,'(a)') 'exit.'
    call exit(0)
  end subroutine lwoniom_write_example

  subroutine write_example_xyz
    integer :: ich,i
    open (newunit=ich,file='example.xyz')
    write (ich,'(2x,i0)') testnat
    write (ich,*)
    do i = 1,testnat
      write (ich,'(a2,1x,3f20.14)') testelem(testat(i)),testxyz(1:3,i)*autoaa
    end do
    close (ich)
  end subroutine write_example_xyz

  subroutine write_example_toml
    integer :: ich,i
    open (newunit=ich,file='example.toml')
    write (ich,'(a)') '# All data for lwONIOM must be contained in a corresponding [lwoniom]-block'
    write (ich,'(a)') '[lwoniom]'
    write (ich,'(a)') '# The systems total number of atoms must be specified'
    write (ich,'(a,i0)') 'natoms = ',testnat
    write (ich,'(a)') '# Then, the XYZ file name can be given. Alternatively, use the CMD argument -i'
    write (ich,'(a)') "xyz = 'example.xyz'"
    write (ich,*)
    write (ich,'(a)') '# Next, fragments must be defined on an by-atom basis'
    write (ich,'(a)') '# An ascending fragment numbering is assumed, i.e., fragment.1 will be the parent system'
    write (ich,'(a)') "fragment.1 = 'all'   # fragment 1 contains all atoms"
    write (ich,'(a)') 'fragment.2 = [1,2]   # fragment 2 contains atoms 1 (Cl) and 2 (Al)'
    write (ich,*)
    write (ich,'(a)') '# Finally, layers are defined on an by-fragment basis'
    write (ich,'(a)') '# As with the fragments, layers are given in ascending order'
    write (ich,'(a)') '# One layer can contain multiple (non-overlapping) fragments in MC-ONIOM, which is not the case here'
    write (ich,'(a)') 'layer.1 = [1]  # layer 1 contains only fragment 1'
    write (ich,'(a)') 'layer.2 = [2]  # layer 2 contains only fragment 2'

    close (ich)
  end subroutine write_example_toml

!========================================================================================!
!========================================================================================!
end module app_example