using OrdinaryDiffEq
using Trixi
using Plots

tau = 0.001
equations = PerturbationMomentSystem1D(0.0, 1.0, tau)

initial_condition = initial_condition_constant
surface_flux = flux_lax_friedrichs

boundary_condition = BoundaryConditionDirichlet(initial_condition)
boundary_conditions = (x_neg=boundary_condition, x_pos=boundary_condition)
volume_flux  = flux_lax_friedrichs
basis = LobattoLegendreBasis(3)
shock_indicator_variable = density_pressure
indicator_sc = IndicatorHennemannGassner(equations, basis,
                                         alpha_max=0.5,
                                         alpha_min=0.001,
                                         alpha_smooth=true,
                                         variable=shock_indicator_variable)
volume_integral = VolumeIntegralShockCapturingHG(indicator_sc;
                                                 volume_flux_dg=volume_flux,
                                                 volume_flux_fv=surface_flux)
solver = DGSEM(basis, surface_flux, volume_integral)
coordinates_min = (-1.0,)
coordinates_max = ( 1.0,)

mesh = TreeMesh(coordinates_min, coordinates_max, initial_refinement_level=9, n_cells_max=10_000, periodicity=false)
semi = SemidiscretizationHyperbolic(mesh, equations, initial_condition, solver, boundary_conditions=boundary_conditions, source_terms=source_terms_convergence_test)

t = 0.2
tspan = (0.0, t)
ode = semidiscretize(semi, tspan)

# At the beginning of the main loop, the SummaryCallback prints a summary of the simulation setup
# and resets the timers
summary_callback = SummaryCallback()

analysis_interval = 100

# The AnalysisCallback allows to analyse the solution in regular intervals and prints the results
analysis_callback = AnalysisCallback(semi, interval=100)

# das hier sorgt für die Ausgabe!
alive_callback = AliveCallback(analysis_interval=analysis_interval)


# The SaveSolutionCallback allows to save the solution to a file in regular intervals
save_solution = SaveSolutionCallback(interval=100,solution_variables=cons2prim)

# The StepsizeCallback handles the re-calculcation of the maximum Δt after each time step
stepsize_callback = StepsizeCallback(cfl=0.9)

save_restart = SaveRestartCallback(interval=100,save_final_restart=true)

# Create a CallbackSet to collect all callbacks such that they can be passed to the ODE solver
callbacks = CallbackSet(summary_callback, analysis_callback, alive_callback, save_solution, stepsize_callback, save_restart)



###############################################################################
# run the simulation

sol = solve(ode, CarpenterKennedy2N54(williamson_condition=false), dt=1.0, save_everystep=false, callback=callbacks);

# Print the timer summary
summary_callback()

pd = PlotData1D(sol; solution_variables=cons2prim)
pd2 = PlotData1D(sol; solution_variables=cons2cons)
# t = 0.4
x =[-0.9966666666666667, -0.99, -0.9833333333333333, -0.9766666666666667, -0.97, -0.9633333333333334, -0.9566666666666667, -0.95, -0.9433333333333334, -0.9366666666666666, -0.9299999999999999, -0.9233333333333333, -0.9166666666666666, -0.91, -0.9033333333333333, -0.8966666666666666, -0.89, -0.8833333333333333, -0.8766666666666667, -0.87, -0.8633333333333333, -0.8566666666666667, -0.85, -0.8433333333333333, -0.8366666666666667, -0.83, -0.8233333333333334, -0.8166666666666667, -0.81, -0.8033333333333333, -0.7966666666666666, -0.79, -0.7833333333333333, -0.7766666666666666, -0.77, -0.7633333333333333, -0.7566666666666666, -0.75, -0.7433333333333333, -0.7366666666666666, -0.73, -0.7233333333333334, -0.7166666666666667, -0.71, -0.7033333333333334, -0.6966666666666667, -0.69, -0.6833333333333333, -0.6766666666666666, -0.6699999999999999, -0.6633333333333333, -0.6566666666666666, -0.6499999999999999, -0.6433333333333333, -0.6366666666666667, -0.6299999999999999, -0.6233333333333333, -0.6166666666666667, -0.61, -0.6033333333333333, -0.5966666666666667, -0.59, -0.5833333333333333, -0.5766666666666667, -0.57, -0.5633333333333332, -0.5566666666666666, -0.55, -0.5433333333333332, -0.5366666666666666, -0.53, -0.5233333333333333, -0.5166666666666666, -0.51, -0.5033333333333333, -0.4966666666666666, -0.49, -0.4833333333333333, -0.4766666666666667, -0.47, -0.46333333333333326, -0.45666666666666667, -0.44999999999999996, -0.44333333333333325, -0.43666666666666665, -0.42999999999999994, -0.42333333333333334, -0.41666666666666663, -0.4099999999999999, -0.4033333333333333, -0.3966666666666666, -0.39, -0.3833333333333333, -0.3766666666666666, -0.37, -0.3633333333333333, -0.3566666666666666, -0.35, -0.34333333333333327, -0.33666666666666667, -0.32999999999999996, -0.32333333333333325, -0.31666666666666665, -0.30999999999999994, -0.30333333333333334, -0.29666666666666663, -0.2899999999999999, -0.2833333333333333, -0.2766666666666666, -0.2699999999999999, -0.2633333333333333, -0.2566666666666666, -0.25, -0.2433333333333333, -0.23666666666666658, -0.22999999999999998, -0.22333333333333327, -0.21666666666666656, -0.20999999999999996, -0.20333333333333325, -0.19666666666666666, -0.18999999999999995, -0.18333333333333324, -0.17666666666666664, -0.16999999999999993, -0.16333333333333333, -0.15666666666666662, -0.1499999999999999, -0.1433333333333333, -0.1366666666666666, -0.1299999999999999, -0.1233333333333333, -0.11666666666666659, -0.10999999999999999, -0.10333333333333328, -0.09666666666666657, -0.08999999999999997, -0.08333333333333326, -0.07666666666666666, -0.06999999999999995, -0.06333333333333324, -0.05666666666666664, -0.04999999999999993, -0.043333333333333224, -0.036666666666666625, -0.029999999999999916, -0.023333333333333317, -0.016666666666666607, -0.009999999999999898, -0.0033333333333332993, 0.0033333333333334103, 0.010000000000000009, 0.01666666666666683, 0.023333333333333428, 0.030000000000000027, 0.036666666666666625, 0.043333333333333446, 0.050000000000000044, 0.05666666666666664, 0.06333333333333346, 0.07000000000000006, 0.07666666666666666, 0.08333333333333348, 0.09000000000000008, 0.09666666666666668, 0.1033333333333335, 0.1100000000000001, 0.1166666666666667, 0.1233333333333333, 0.13000000000000012, 0.13666666666666671, 0.1433333333333333, 0.15000000000000013, 0.15666666666666673, 0.16333333333333333, 0.17000000000000015, 0.17666666666666675, 0.18333333333333335, 0.19000000000000017, 0.19666666666666677, 0.20333333333333337, 0.2100000000000002, 0.21666666666666679, 0.22333333333333338, 0.22999999999999998, 0.2366666666666668, 0.2433333333333334, 0.25, 0.2566666666666668, 0.2633333333333334, 0.27, 0.27666666666666684, 0.28333333333333344, 0.29000000000000004, 0.29666666666666686, 0.30333333333333345, 0.31000000000000005, 0.31666666666666665, 0.3233333333333335, 0.33000000000000007, 0.33666666666666667, 0.3433333333333335, 0.3500000000000001, 0.3566666666666667, 0.3633333333333335, 0.3700000000000001, 0.3766666666666667, 0.3833333333333335, 0.3900000000000001, 0.3966666666666667, 0.4033333333333333, 0.41000000000000014, 0.41666666666666674, 0.42333333333333334, 0.43000000000000016, 0.43666666666666676, 0.44333333333333336, 0.4500000000000002, 0.4566666666666668, 0.4633333333333334, 0.4700000000000002, 0.4766666666666668, 0.4833333333333334, 0.49, 0.4966666666666668, 0.5033333333333334, 0.51, 0.5166666666666668, 0.5233333333333334, 0.53, 0.5366666666666668, 0.5433333333333334, 0.55, 0.5566666666666669, 0.5633333333333335, 0.5700000000000001, 0.5766666666666667, 0.5833333333333335, 0.5900000000000001, 0.5966666666666667, 0.6033333333333335, 0.6100000000000001, 0.6166666666666667, 0.6233333333333335, 0.6300000000000001, 0.6366666666666667, 0.6433333333333335, 0.6500000000000001, 0.6566666666666667, 0.6633333333333333, 0.6700000000000002, 0.6766666666666667, 0.6833333333333333, 0.6900000000000002, 0.6966666666666668, 0.7033333333333334, 0.7100000000000002, 0.7166666666666668, 0.7233333333333334, 0.7300000000000002, 0.7366666666666668, 0.7433333333333334, 0.7500000000000002, 0.7566666666666668, 0.7633333333333334, 0.77, 0.7766666666666668, 0.7833333333333334, 0.79, 0.7966666666666669, 0.8033333333333335, 0.81, 0.8166666666666669, 0.8233333333333335, 0.8300000000000001, 0.8366666666666669, 0.8433333333333335, 0.8500000000000001, 0.8566666666666667, 0.8633333333333335, 0.8700000000000001, 0.8766666666666667, 0.8833333333333335, 0.8900000000000001, 0.8966666666666667, 0.9033333333333335, 0.9100000000000001, 0.9166666666666667, 0.9233333333333336, 0.9300000000000002, 0.9366666666666668, 0.9433333333333334, 0.9500000000000002, 0.9566666666666668, 0.9633333333333334, 0.9700000000000002, 0.9766666666666668, 0.9833333333333334, 0.9900000000000002, 0.9966666666666668];
rho = [2.9999999999971596, 2.999999999987033, 2.999999999957886, 2.999999999859364, 2.9999999995330704, 2.999999998453627, 2.9999999948917564, 2.9999999831644186, 2.9999999446302397, 2.9999998182437255, 2.9999994044007163, 2.9999980513449267, 2.999993633604031, 2.999979227453837, 2.9999323024301914, 2.9997796603281515, 2.9992838931655554, 2.9976749132233196, 2.99263843731171, 2.98217310401934, 2.966907006099741, 2.948172611598559, 2.927753160656146, 2.907882530935729, 2.8906130968009465, 2.8781072019896135, 2.8710615380089504, 2.8680337964429707, 2.8671010780754287, 2.8667567448725535, 2.8666491485920194, 2.8666298440482527, 2.8666292268447937, 2.8666437325025758, 2.8666964117511333, 2.866763583773931, 2.8668186021370596, 2.866856671853383, 2.8668903843745257, 2.866935276679551, 2.8669983337998635, 2.867065735346075, 2.8671148021902564, 2.867136348535113, 2.8671384120338677, 2.8671327378542975, 2.8670922372487233, 2.867025055835664, 2.866970284272035, 2.866947169595476, 2.8669453091086385, 2.8669581299521987, 2.867028747681405, 2.8671276726802883, 2.8672020975146264, 2.8672292048781234, 2.8672307191099287, 2.867212980004488, 2.867129902356469, 2.8670359359083557, 2.866974931177554, 2.866966696416921, 2.8669720901735376, 2.8670266045851216, 2.8671466425696543, 2.8672531674273083, 2.867288609932804, 2.8672848760402747, 2.867258304739503, 2.8671373910045395, 2.8669778958202814, 2.866873795862651, 2.8668710859512503, 2.8668937362101774, 2.8669988449296055, 2.867217877510279, 2.8674546083184116, 2.8675570991133053, 2.8675468916178852, 2.8674856309839476, 2.867257005717766, 2.866894064362711, 2.8665403179548603, 2.866368936912509, 2.8663529390051963, 2.8664442560729264, 2.866783866691608, 2.8673217381786666, 2.8679180422964707, 2.8684609585327925, 2.868782375151878, 2.8688957614782327, 2.868687436986887, 2.868086150289442, 2.866975912318968, 2.8646798511190634, 2.858471803231329, 2.8388046936081563, 2.794286564800019, 2.7253147011654177, 2.6384250122390625, 2.5426674762085515, 2.449262959529336, 2.369463226806091, 2.3122892521514444, 2.283050333978335, 2.270877736449921, 2.2661398792013117, 2.2646365007195577, 2.264415072441482, 2.264172966387209, 2.2641556432663523, 2.2641526812417814, 2.264398548082754, 2.2649516200713413, 2.265520717755241, 2.265920517352683, 2.2661124220075326, 2.2662022870874168, 2.2664171997500953, 2.2666771940569213, 2.266955430531363, 2.267195065288975, 2.2673440197309276, 2.2673518222571007, 2.2673308247750157, 2.2673009131470803, 2.2672482338816566, 2.2670670984298362, 2.2666922156326996, 2.2661636492018538, 2.2657337051235316, 2.265553842004543, 2.2655307159180857, 2.2655938039722843, 2.26579494688254, 2.266395591833907, 2.267530315673026, 2.2687399432098094, 2.269474559281484, 2.2696858933696826, 2.2696847800188578, 2.269487150649506, 2.2686883572651593, 2.2660458086593165, 2.2582631088265734, 2.236502991095803, 2.190955403621243, 2.1241356754780165, 2.043066614858798, 1.9569333851373547, 1.8758643245184075, 1.8090445963753141, 1.763497008900806, 1.7417368911722835, 1.733954191340176, 1.7313116427352884, 1.7305128493515511, 1.730315219981628, 1.7303141066302856, 1.7305254407194954, 1.7312600567917453, 1.7324696843283347, 1.73360440816738, 1.7342050531170665, 1.7344061960270123, 1.7344692840812908, 1.7344461579948052, 1.7342662948774115, 1.7338363508008645, 1.7333077843696465, 1.7329329015718171, 1.7327517661199126, 1.7326990868549703, 1.7326691752271244, 1.732648177744926, 1.7326559802694121, 1.7328049347085135, 1.7330445694657173, 1.733322805942058, 1.733582800250575, 1.7337977129176105, 1.733887577998848, 1.7340794826498784, 1.7344792822432307, 1.735048379925165, 1.735601451915893, 1.7358473187607164, 1.7358443567385913, 1.735827033617042, 1.7355849275537545, 1.7353634992678093, 1.7338601207925701, 1.7291222635625438, 1.716949666016267, 1.6877107478488225, 1.630536773193283, 1.5507370404733194, 1.4573325237925119, 1.3615749877623482, 1.2746852988347757, 1.205713435198065, 1.1611953063895466, 1.1415281967726116, 1.1353201488838707, 1.133024087681783, 1.1319138497096726, 1.1313125630112415, 1.1311042385200525, 1.131217624846748, 1.1315390414667301, 1.1320819577037695, 1.1326782618214672, 1.1332161333080584, 1.1335557439268236, 1.1336470609940874, 1.1336310630879893, 1.133459682047254, 1.1331059356399402, 1.1327429942848692, 1.1325143690162742, 1.1324531083818512, 1.1324429008827455, 1.1325453916767025, 1.1327821224846675, 1.1330011550691792, 1.1331062637904448, 1.1331289140523757, 1.1331262041387824, 1.1330221041773214, 1.1328626089901659, 1.1327416952569926, 1.1327151239571922, 1.1327113900720063, 1.1327468325789778, 1.1328533574336883, 1.1329733954138916, 1.1330279098249783, 1.1330333035797402, 1.1330250688226426, 1.1329640640906482, 1.132870097642622, 1.1327870199911532, 1.1327692808848875, 1.1327707951176493, 1.132797902485734, 1.1328723273220924, 1.1329712523217748, 1.1330418700479132, 1.1330546908913088, 1.1330528304037712, 1.133029715730126, 1.132974944167338, 1.1329077627522597, 1.1328672621431424, 1.1328615879630164, 1.1328636514612385, 1.1328851978086538, 1.1329342646550424, 1.1330016662035525, 1.1330647233263678, 1.1331096156320772, 1.1331433281508334, 1.1331813978585408, 1.133236416218356, 1.1333035882425722, 1.1333562675031317, 1.133370773162378, 1.1333701559644191, 1.1333508513964206, 1.1332432551191398, 1.132898921925896, 1.131966203563651, 1.1289384619946887, 1.1218927980085103, 1.109386903193744, 1.0921174690619964, 1.0722468393446496, 1.051827388407045, 1.0330929938990197, 1.0178268959876138, 1.007361562677124, 1.0023250867817173, 1.000716106832348, 1.0002203396724485, 1.000067697569734, 1.0000207725461103, 1.0000063663960157, 1.0000019486550318, 1.0000005955992908, 1.0000001817562627, 1.0000000553697603, 1.0000000168355723, 1.0000000051082387, 1.0000000015463644, 1.0000000004669227, 1.0000000001406277, 1.0000000000421079, 1.0000000000129612, 1.0000000000028355];
v = [2.0158084431138644*10^-12, 9.207532775441496*10^-12, 2.990648522248397*10^-11, 9.987394852179215*10^-11, 3.315970458440531*10^-10, 1.0981852568670473*10^-9, 3.62772057744507*10^-9, 1.1956118737480355*10^-8, 3.932194257034089*10^-8, 1.2907783806371224*10^-7, 4.2297679508241717*10^-7, 1.3838771159806923*10^-6, 4.521232554733502*10^-6, 0.000014752139881268724, 0.00004807786270662424, 0.00015649010414744765, 0.0005086789582158509, 0.0016524874440453307, 0.005240820643439441, 0.012735805450930327, 0.023763787482101766, 0.03745328781385868, 0.05257362417230073, 0.0674914334798327, 0.08062300715077811, 0.09023079275424906, 0.09568054101230232, 0.09803055142197506, 0.09875550362187814, 0.09902327106279234, 0.09910703923897808, 0.09912218395620295, 0.09912270992324586, 0.0991114742343799, 0.0990705561123924, 0.09901833433946876, 0.09897552168911403, 0.09894585395702608, 0.09891954170378016, 0.09888451178446926, 0.09883525458346866, 0.09878257924669957, 0.09874415171759571, 0.0987276241613436, 0.09872617837814596, 0.09873086874166755, 0.09876291392277255, 0.09881584516546149, 0.09885917093965996, 0.09887694484613309, 0.09887816768834591, 0.09886779234924602, 0.09881204051179113, 0.09873396478606335, 0.09867464075277574, 0.09865310976569994, 0.09865192870306941, 0.09866605281257212, 0.09873156081375103, 0.09880638845660963, 0.09885649920396264, 0.09886514978965538, 0.09886170698420935, 0.09882068857148749, 0.09872700750291609, 0.0986414685141015, 0.09860470778526836, 0.09860410991202397, 0.09862069403828393, 0.09870848717532184, 0.09882945795742948, 0.098915967685019, 0.09893288091984545, 0.09892187396278371, 0.09886425369300755, 0.09871098345475003, 0.09854994480103785, 0.09846182981909082, 0.09845482462415135, 0.09848651802696921, 0.0986239381974604, 0.09884912030089171, 0.09905224638955659, 0.09914901985741367, 0.0991572295345744, 0.09909469620232994, 0.09887994589189866, 0.09857571816553926, 0.09828559767371416, 0.09808285390609836, 0.09802064126462647, 0.09806390871000259, 0.09823081988485517, 0.098638396808528, 0.09919169237623565, 0.09980694289888352, 0.10139062224008169, 0.10751596320085288, 0.12248529902894828, 0.14615398497884982, 0.17771871938174505, 0.21505028673870602, 0.25425592527855384, 0.29024491130366975, 0.31748369493659023, 0.33248030849877, 0.33818883072653794, 0.3400098109655885, 0.3411482716758572, 0.3417555414680841, 0.3418897395468298, 0.34173538311365337, 0.34148580746771284, 0.3412221093144884, 0.34098113508541344, 0.34083037657408605, 0.3406977967355811, 0.34048531989014164, 0.3403072718480582, 0.34021603722821675, 0.3401976996580129, 0.34019903057243167, 0.34019184116324525, 0.3401080918494081, 0.3400074162523882, 0.33995308005953156, 0.3399737315970089, 0.3400518932783367, 0.3401613089977195, 0.3402902086238418, 0.34038066843592935, 0.3404298143868475, 0.3404418001099148, 0.3404445677212481, 0.3404141238932229, 0.34036246930482067, 0.34027191119494166, 0.34010707962012604, 0.3399457286304841, 0.33985970274806515, 0.3398393957528415, 0.33982275047195803, 0.33983940557959735, 0.3399131795690188, 0.3402692678287254, 0.34142666127279253, 0.34477450084789385, 0.3520125367900526, 0.36315761348384296, 0.37760851717959826, 0.3942287258187992, 0.4112216606167194, 0.42632656550600434, 0.4372500766975734, 0.4426795099947918, 0.44468634294545717, 0.44541782885210823, 0.44568358133846026, 0.4457514536918534, 0.4457737930909296, 0.445704484290927, 0.44548388327611366, 0.4451466715885607, 0.4448481763880249, 0.4446945657726627, 0.4446709955508917, 0.4446822047056977, 0.44468905803880676, 0.4447548262689432, 0.44488529574435115, 0.4450064633051437, 0.4450077155829679, 0.4449466274047677, 0.44486821626749046, 0.4448547411226438, 0.44493535659645206, 0.44506356535843916, 0.4451056481278198, 0.4450064663148972, 0.4448786830967351, 0.44478491497054284, 0.4448068606975094, 0.44499886999359467, 0.4451895864995906, 0.4451816099980076, 0.44512059909658, 0.4451844909737679, 0.4454170584141754, 0.4457439361804323, 0.4459531340210241, 0.44588794642831653, 0.4451959653067702, 0.4443898228649741, 0.4441475901203124, 0.44210339672767374, 0.43497627568192876, 0.4217780643793024, 0.40157654313872754, 0.3752070038352845, 0.34437876619228774, 0.3124815232905431, 0.2838643208776972, 0.26284710185164695, 0.25388968543302887, 0.2518363992853223, 0.2509921861675326, 0.24993370285639288, 0.24908546774749926, 0.2487260877193028, 0.24858160081900701, 0.2486408571152866, 0.2489881911364209, 0.2495397934980624, 0.25014419164462354, 0.25058266808540114, 0.2507126124866561, 0.2506967918358935, 0.25050494724273037, 0.25010014274720477, 0.249641956869172, 0.2493643197961499, 0.24930288435767045, 0.24932375740937598, 0.24951537965617784, 0.24985024998359207, 0.25017070801408026, 0.25028411712171555, 0.25030498492827424, 0.2502630286994595, 0.2500762080295029, 0.24981916796933212, 0.24963414447454152, 0.24960033382748623, 0.2496029950671824, 0.2496851501922758, 0.2498688874658905, 0.25006901694798395, 0.25015602197929543, 0.25016307198990634, 0.2501437195021349, 0.2500357031485098, 0.24987526011112626, 0.24973501843441318, 0.24970472386548326, 0.24970724766001717, 0.2497534081857207, 0.24988065807209392, 0.25004779262906224, 0.2501671196434387, 0.2501894227230246, 0.2501869017585666, 0.25014904858603665, 0.2500562836423699, 0.2499430166725619, 0.24987438111151675, 0.249864256422035, 0.24986728056617896, 0.2499024787039258, 0.24998462579786035, 0.25009716989505887, 0.2502023841219658, 0.2502771832602623, 0.2503333642935539, 0.25039668607602905, 0.2504880278875563, 0.25059940750038184, 0.25068665040843735, 0.2507105918214065, 0.2507094519927837, 0.2506771043811129, 0.25049840707303056, 0.24992698414333814, 0.24837749895675162, 0.24333011096104606, 0.23147835063103053, 0.2100709136861519, 0.1797033432515434, 0.14355145540589986, 0.10497802068735211, 0.06824646763904217, 0.03731516293005586, 0.01556926718719355, 0.0049421292648777265, 0.001524580847686992, 0.0004693324176037691, 0.00014422057017438064, 0.00004425519403873766, 0.000013563582408070524, 4.151620616227469*10^-6, 1.2689293388954171*10^-6, 3.872334331382588*10^-7, 1.1796581445615875*10^-7, 3.58683517091927*10^-8, 1.0883157467740594*10^-8, 3.294546808585163*10^-9, 9.947857873552446*10^-10, 2.9961078181934263*10^-10, 8.971335691411621*10^-11, 2.761585640480297*10^-11, 6.042414052703605*10^-12];
theta = [0.9999999999985604, 0.9999999999934237, 0.9999999999786398, 0.9999999999286657, 0.9999999997631587, 0.9999999992156261, 0.9999999974089187, 0.9999999914603964, 0.9999999719144901, 0.999999907806748, 0.9999996978907547, 0.9999990115716032, 0.9999967707240061, 0.9999894632807286, 0.9999656598408934, 0.9998882195141148, 0.9996365921823837, 0.9988188083502382, 0.9962476129852594, 0.9908494436615837, 0.9828385868136218, 0.9727815819304663, 0.9615282077211619, 0.9502761979963682, 0.9402486881777106, 0.9328391830044899, 0.9286090118566832, 0.9267789214571222, 0.9262136026380597, 0.9260046949885761, 0.9259392504409435, 0.9259273168582778, 0.9259268620176127, 0.925935586333658, 0.9259674655024491, 0.9260081941562629, 0.9260416210264746, 0.9260648245161313, 0.9260854414453444, 0.9261128821401234, 0.9261515169239073, 0.9261928560641349, 0.926223087627799, 0.9262357685131615, 0.9262367501563431, 0.926232832659591, 0.9262073167534786, 0.9261653980948591, 0.9261309325594653, 0.926117213476613, 0.9261164416221959, 0.9261249193841997, 0.9261691712647724, 0.9262312526235895, 0.9262788200804484, 0.9262960744640669, 0.9262969664558687, 0.9262856437179352, 0.9262336831797262, 0.9261737049902751, 0.9261321164593608, 0.9261232827441553, 0.9261252204279572, 0.926156221047698, 0.9262297189661529, 0.9262984711635174, 0.9263353489171021, 0.9263397256373369, 0.9263303887217229, 0.9262677007001995, 0.9261764972893276, 0.926105505223763, 0.9260761193538063, 0.926074353767592, 0.9261075393343414, 0.9262052225105265, 0.92631148782949, 0.9263807687168477, 0.9264044611430916, 0.9263965346368055, 0.9263250680858615, 0.9262037447798701, 0.9261074382506974, 0.9260669149159245, 0.9260674027642136, 0.9261104290245651, 0.9262409479274524, 0.926357198157115, 0.9264122558862511, 0.9263755902005061, 0.9262760327463372, 0.9261391180682065, 0.9259230502767168, 0.9257680964829554, 0.9257449732559804, 0.926239016443396, 0.9279789746889879, 0.9317775960717998, 0.9390905375043959, 0.9515435342653706, 0.9677006126154484, 0.9858242729042235, 1.0038733341419244, 1.0195794755291931, 1.0309630303340591, 1.036919150453007, 1.0390538873829942, 1.0397115700679609, 1.0401709645827577, 1.0404689245755991, 1.0405186523311667, 1.0404650680461283, 1.0403389801749292, 1.0401852568723402, 1.0400905711353245, 1.0400257013930012, 1.0399599365943284, 1.039846158488334, 1.0397482400390408, 1.0397072010725836, 1.0397270490294723, 1.03979588877381, 1.03985282262186, 1.0398992370584101, 1.039921117252011, 1.0399306915313167, 1.0399700348644585, 1.0399866661453745, 1.0399606151651928, 1.039875848161901, 1.0397145025834107, 1.03955711249424, 1.0394777438305647, 1.0394609869675453, 1.0394752014585749, 1.0395242230591164, 1.0396911477349482, 1.0400304770677378, 1.040398366603338, 1.040645249597374, 1.0407299244858883, 1.0407433294470696, 1.0406884320730976, 1.040452398482782, 1.0396653825572355, 1.0373512920462387, 1.0308270747033477, 1.0167545770322486, 0.9949711636605689, 0.966515218568059, 0.9335317893838722, 0.8995472277082256, 0.8690994539345964, 0.8469243709094398, 0.8358693004509763, 0.8318096461521952, 0.8303915566273317, 0.8299407802581606, 0.8298326466611975, 0.8298385270647493, 0.8299865013727168, 0.8304276445806281, 0.8310887040348953, 0.8316907672563877, 0.8319901261469368, 0.8320667390381957, 0.8320795096597348, 0.8320540090491824, 0.8319156899248552, 0.8316411458976287, 0.8313659388310325, 0.8312537442527923, 0.8312512917381388, 0.831313983579581, 0.8313724601872969, 0.8313427678539208, 0.8313043136358956, 0.8313432765192508, 0.8314700791684918, 0.8316274376723567, 0.8317028574152274, 0.8316380494371194, 0.8314096858284299, 0.8311614174574834, 0.831080363354603, 0.8310286270522056, 0.8308725629089169, 0.830550525916504, 0.8302145135179526, 0.8300347331199593, 0.830133341874073, 0.8308847618184214, 0.8319042568675213, 0.8328857661796605, 0.8367594290620098, 0.8484776002881951, 0.87144227600928, 0.9060934054138543, 0.9509100792407682, 1.0026558266504035, 1.0558296430872274, 1.1027106006032823, 1.1343354850358716, 1.150278687523788, 1.1565974514550754, 1.158595263220131, 1.1590516346176518, 1.1590009645698451, 1.1585864165806774, 1.158245277826308, 1.1579019359824088, 1.1575984379185056, 1.157467064380875, 1.1574918173834856, 1.1576382865212822, 1.1576948412258754, 1.1577036092791162, 1.1576870413740774, 1.157622632610882, 1.1575091628707983, 1.1574474632988678, 1.1574563176350854, 1.1575139179102065, 1.1576203148984572, 1.1577547437165439, 1.1578748130876284, 1.157907114269571, 1.157892342240737, 1.1578284135832397, 1.157717272255174, 1.1575845131663793, 1.1574978202961284, 1.157488405689977, 1.1574997621229823, 1.1575655363023238, 1.157671307633697, 1.1577819534875684, 1.157827058969954, 1.1578290252634253, 1.157813154687292, 1.157747699200202, 1.1576557752372074, 1.1575767830105645, 1.1575594148360198, 1.1575607232795466, 1.1575871720251725, 1.1576601368814408, 1.1577546337589, 1.1578218833601306, 1.1578350632949275, 1.1578341441131896, 1.1578137937540847, 1.157761320228694, 1.157697711725199, 1.1576588126963159, 1.1576526129784517, 1.157653899296291, 1.1576726398847235, 1.1577182697187975, 1.157780516417501, 1.1578386419979925, 1.1578798694363188, 1.1579108444120003, 1.157945644693477, 1.1579956916822158, 1.1580565930981679, 1.1581042028840849, 1.1581171923872222, 1.1581164564587876, 1.158098465509001, 1.1580004793227339, 1.1576876533743152, 1.1568383841731753, 1.154060931749529, 1.1474716638998188, 1.1353322804052695, 1.1175878333018008, 1.095661944610056, 1.0713065972817897, 1.04719219331211, 1.0261880400040928, 1.0110394754883558, 1.0035217511314778, 1.0010881499496371, 1.0003351450553153, 1.0001030019344017, 1.0000316084118048, 1.000009687664062, 1.000002965269771, 1.0000009063263193, 1.000000276579613, 1.0000000842565249, 1.0000000256188037, 1.0000000077732438, 1.0000000023531153, 1.0000000007105216, 1.000000000213996, 1.0000000000640779, 1.0000000000197247, 1.0000000000043159];


#rho plot
plot(pd.x, pd.data[:,1], xlims = (-2.0, 2.0), label = "DGSEM ohne amr", title ="rho, t = "*string(t)*", tau = "*string(tau)*", B ="*string(1/tau), seriestype = :scatter, markersize=1)
#plot!(x,rho, label = "FVV")



#vx plot
# plot(pd.x, pd.data[:,2], xlims = (-1.1,1.1), label = "DGSEM", title ="rho, t = "*string(t)*"", tau = "*string(tau)*", B ="*string(1/tau), seriestype = :scatter, markersize=1)
#plot!(x,v, label = "FVV")


# # #theta plot
#plot(pd.x, pd.data[:,3], xlims = (-1.1,1.1), ylims = (0.5, 1.5), label = "DGSEM", title ="theta, t = 0.4", seriestype = :scatter, markersize=1)
#plot!(x,theta, label = "FVV")



#savefig("rho.png")

# bildx =[-1.70073327222731,	-1.63473877176902,	-1.55499541704858,	-1.480751604033,	-1.40650779101742,	-1.33776351970669,	-1.26076993583868,	-1.17277726856095,	-1.07103574702108,	-0.972043996333639,	-0.867552703941338,	-0.774060494958754,	-0.67781851512374,	-0.578826764436297,	-0.507332722273144,	-0.441338221814849,	-0.378093492208983,	-0.325847846012832,	-0.270852428964253,	-0.210357470210816,	-0.152612282309808,	-0.0976168652612281,	-0.037121906507791,	0.00412465627864345,	0.0371219065077908,	0.10724106324473,	0.185609532538955,	0.262603116406966,	0.336846929422548,	0.400091659028414,	0.457836846929423,	0.488084326306141,	0.504582951420715,	0.559578368469294,	0.644821264894593,	0.768560953253895,	0.897800183318057,	1.04353803849679,	1.06278643446379,	1.07378551787351,	1.0820348304308,	1.16452795600367,	1.24977085242896,	1.34601283226398,	1.41200733272227,	1.50274977085243,	1.57149404216315,	1.64573785517874,	1.74472960586618,	1.76397800183318]
# bildy =[5.00679873273824,	5.01254049650017,	4.99420318577056,	4.98189390123396,	4.94556059267333,	4.8912312862817,	4.84287495561556,	4.7824625725634,	4.69196510901369,	4.60748466156349,	4.52898819993229,	4.44452977266909,	4.32402429927636,	4.22152583380814,	4.13715548729298,	4.07082517898283,	3.96246383872415,	3.7760684607614,	3.61368609672918,	3.42725768848592,	3.21081026030614,	3.00037984822586,	2.78992741595858,	2.60357607836985,	2.48332383712769,	2.42298301968329,	2.34759415694338,	2.28422281126222,	2.22686848168058,	2.19057921349397,	2.16932698050938,	2.03707373734871,	1.83280347258349,	1.82057125870141,	1.78419390976678,	1.71763238949307,	1.65705485503836,	1.57839324200461,	1.38612397915789,	1.19989375259769,	1.1157766382331,	1.10343432341599,	1.10009000751439,	1.09670165123877,	1.09043140298869,	1.09006806990308,	1.0867898145625,	1.08348953903491,	1.07708716966279,	1.02896205096022]
# plot!(bildx, bildy, label ="Paper")

