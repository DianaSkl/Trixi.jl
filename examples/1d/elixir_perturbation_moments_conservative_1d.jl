using OrdinaryDiffEq
using Trixi
using Plots

equations = PerturbationMomentSystem1D(0.0, 1.0)

initial_condition = initial_condition_constant
# initial_condition = initial_condition_convergence_test
boundary_condition = BoundaryConditionDirichlet(initial_condition)
boundary_conditions = (x_neg=boundary_condition,
                       x_pos=boundary_condition)

solver = DGSEM(polydeg=3, surface_flux=flux_lax_friedrichs)

coordinates_min = (-1,)
coordinates_max = ( 1,)

mesh = TreeMesh(coordinates_min, coordinates_max, initial_refinement_level=7, n_cells_max=10_000, periodicity=false)
#semi = SemidiscretizationHyperbolic(mesh, equations, initial_condition, solver)
semi = SemidiscretizationHyperbolic(mesh, equations, initial_condition, solver, boundary_conditions=boundary_conditions)

tspan = (0.0, 0.4)
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


# t = 0.05
x =[-0.9966666666666667, -0.99, -0.9833333333333333, -0.9766666666666667, -0.97, -0.9633333333333334, -0.9566666666666667, -0.95, -0.9433333333333334, -0.9366666666666666, -0.9299999999999999, -0.9233333333333333, -0.9166666666666666, -0.91, -0.9033333333333333, -0.8966666666666666, -0.89, -0.8833333333333333, -0.8766666666666667, -0.87, -0.8633333333333333, -0.8566666666666667, -0.85, -0.8433333333333333, -0.8366666666666667, -0.83, -0.8233333333333334, -0.8166666666666667, -0.81, -0.8033333333333333, -0.7966666666666666, -0.79, -0.7833333333333333, -0.7766666666666666, -0.77, -0.7633333333333333, -0.7566666666666666, -0.75, -0.7433333333333333, -0.7366666666666666, -0.73, -0.7233333333333334, -0.7166666666666667, -0.71, -0.7033333333333334, -0.6966666666666667, -0.69, -0.6833333333333333, -0.6766666666666666, -0.6699999999999999, -0.6633333333333333, -0.6566666666666666, -0.6499999999999999, -0.6433333333333333, -0.6366666666666667, -0.6299999999999999, -0.6233333333333333, -0.6166666666666667, -0.61, -0.6033333333333333, -0.5966666666666667, -0.59, -0.5833333333333333, -0.5766666666666667, -0.57, -0.5633333333333332, -0.5566666666666666, -0.55, -0.5433333333333332, -0.5366666666666666, -0.53, -0.5233333333333333, -0.5166666666666666, -0.51, -0.5033333333333333, -0.4966666666666666, -0.49, -0.4833333333333333, -0.4766666666666667, -0.47, -0.46333333333333326, -0.45666666666666667, -0.44999999999999996, -0.44333333333333325, -0.43666666666666665, -0.42999999999999994, -0.42333333333333334, -0.41666666666666663, -0.4099999999999999, -0.4033333333333333, -0.3966666666666666, -0.39, -0.3833333333333333, -0.3766666666666666, -0.37, -0.3633333333333333, -0.3566666666666666, -0.35, -0.34333333333333327, -0.33666666666666667, -0.32999999999999996, -0.32333333333333325, -0.31666666666666665, -0.30999999999999994, -0.30333333333333334, -0.29666666666666663, -0.2899999999999999, -0.2833333333333333, -0.2766666666666666, -0.2699999999999999, -0.2633333333333333, -0.2566666666666666, -0.25, -0.2433333333333333, -0.23666666666666658, -0.22999999999999998, -0.22333333333333327, -0.21666666666666656, -0.20999999999999996, -0.20333333333333325, -0.19666666666666666, -0.18999999999999995, -0.18333333333333324, -0.17666666666666664, -0.16999999999999993, -0.16333333333333333, -0.15666666666666662, -0.1499999999999999, -0.1433333333333333, -0.1366666666666666, -0.1299999999999999, -0.1233333333333333, -0.11666666666666659, -0.10999999999999999, -0.10333333333333328, -0.09666666666666657, -0.08999999999999997, -0.08333333333333326, -0.07666666666666666, -0.06999999999999995, -0.06333333333333324, -0.05666666666666664, -0.04999999999999993, -0.043333333333333224, -0.036666666666666625, -0.029999999999999916, -0.023333333333333317, -0.016666666666666607, -0.009999999999999898, -0.0033333333333332993, 0.0033333333333334103, 0.010000000000000009, 0.01666666666666683, 0.023333333333333428, 0.030000000000000027, 0.036666666666666625, 0.043333333333333446, 0.050000000000000044, 0.05666666666666664, 0.06333333333333346, 0.07000000000000006, 0.07666666666666666, 0.08333333333333348, 0.09000000000000008, 0.09666666666666668, 0.1033333333333335, 0.1100000000000001, 0.1166666666666667, 0.1233333333333333, 0.13000000000000012, 0.13666666666666671, 0.1433333333333333, 0.15000000000000013, 0.15666666666666673, 0.16333333333333333, 0.17000000000000015, 0.17666666666666675, 0.18333333333333335, 0.19000000000000017, 0.19666666666666677, 0.20333333333333337, 0.2100000000000002, 0.21666666666666679, 0.22333333333333338, 0.22999999999999998, 0.2366666666666668, 0.2433333333333334, 0.25, 0.2566666666666668, 0.2633333333333334, 0.27, 0.27666666666666684, 0.28333333333333344, 0.29000000000000004, 0.29666666666666686, 0.30333333333333345, 0.31000000000000005, 0.31666666666666665, 0.3233333333333335, 0.33000000000000007, 0.33666666666666667, 0.3433333333333335, 0.3500000000000001, 0.3566666666666667, 0.3633333333333335, 0.3700000000000001, 0.3766666666666667, 0.3833333333333335, 0.3900000000000001, 0.3966666666666667, 0.4033333333333333, 0.41000000000000014, 0.41666666666666674, 0.42333333333333334, 0.43000000000000016, 0.43666666666666676, 0.44333333333333336, 0.4500000000000002, 0.4566666666666668, 0.4633333333333334, 0.4700000000000002, 0.4766666666666668, 0.4833333333333334, 0.49, 0.4966666666666668, 0.5033333333333334, 0.51, 0.5166666666666668, 0.5233333333333334, 0.53, 0.5366666666666668, 0.5433333333333334, 0.55, 0.5566666666666669, 0.5633333333333335, 0.5700000000000001, 0.5766666666666667, 0.5833333333333335, 0.5900000000000001, 0.5966666666666667, 0.6033333333333335, 0.6100000000000001, 0.6166666666666667, 0.6233333333333335, 0.6300000000000001, 0.6366666666666667, 0.6433333333333335, 0.6500000000000001, 0.6566666666666667, 0.6633333333333333, 0.6700000000000002, 0.6766666666666667, 0.6833333333333333, 0.6900000000000002, 0.6966666666666668, 0.7033333333333334, 0.7100000000000002, 0.7166666666666668, 0.7233333333333334, 0.7300000000000002, 0.7366666666666668, 0.7433333333333334, 0.7500000000000002, 0.7566666666666668, 0.7633333333333334, 0.77, 0.7766666666666668, 0.7833333333333334, 0.79, 0.7966666666666669, 0.8033333333333335, 0.81, 0.8166666666666669, 0.8233333333333335, 0.8300000000000001, 0.8366666666666669, 0.8433333333333335, 0.8500000000000001, 0.8566666666666667, 0.8633333333333335, 0.8700000000000001, 0.8766666666666667, 0.8833333333333335, 0.8900000000000001, 0.8966666666666667, 0.9033333333333335, 0.9100000000000001, 0.9166666666666667, 0.9233333333333336, 0.9300000000000002, 0.9366666666666668, 0.9433333333333334, 0.9500000000000002, 0.9566666666666668, 0.9633333333333334, 0.9700000000000002, 0.9766666666666668, 0.9833333333333334, 0.9900000000000002, 0.9966666666666668];

w0 = [2.9999999999971596, 2.999999999987033, 2.999999999957886, 2.999999999859364, 2.9999999995330704, 2.999999998453627, 2.9999999948917564, 2.9999999831644186, 2.9999999446302397, 2.9999998182437255, 2.9999994044007163, 2.9999980513449267, 2.999993633604031, 2.999979227453837, 2.9999323024301914, 2.9997796603281515, 2.9992838931655554, 2.9976749132233196, 2.99263843731171, 2.98217310401934, 2.966907006099741, 2.948172611598559, 2.927753160656146, 2.907882530935729, 2.8906130968009465, 2.8781072019896135, 2.8710615380089504, 2.8680337964429707, 2.8671010780754287, 2.8667567448725535, 2.8666491485920194, 2.8666298440482527, 2.8666292268447937, 2.8666437325025758, 2.8666964117511333, 2.866763583773931, 2.8668186021370596, 2.866856671853383, 2.8668903843745257, 2.866935276679551, 2.8669983337998635, 2.867065735346075, 2.8671148021902564, 2.867136348535113, 2.8671384120338677, 2.8671327378542975, 2.8670922372487233, 2.867025055835664, 2.866970284272035, 2.866947169595476, 2.8669453091086385, 2.8669581299521987, 2.867028747681405, 2.8671276726802883, 2.8672020975146264, 2.8672292048781234, 2.8672307191099287, 2.867212980004488, 2.867129902356469, 2.8670359359083557, 2.866974931177554, 2.866966696416921, 2.8669720901735376, 2.8670266045851216, 2.8671466425696543, 2.8672531674273083, 2.867288609932804, 2.8672848760402747, 2.867258304739503, 2.8671373910045395, 2.8669778958202814, 2.866873795862651, 2.8668710859512503, 2.8668937362101774, 2.8669988449296055, 2.867217877510279, 2.8674546083184116, 2.8675570991133053, 2.8675468916178852, 2.8674856309839476, 2.867257005717766, 2.866894064362711, 2.8665403179548603, 2.866368936912509, 2.8663529390051963, 2.8664442560729264, 2.866783866691608, 2.8673217381786666, 2.8679180422964707, 2.8684609585327925, 2.868782375151878, 2.8688957614782327, 2.868687436986887, 2.868086150289442, 2.866975912318968, 2.8646798511190634, 2.858471803231329, 2.8388046936081563, 2.794286564800019, 2.7253147011654177, 2.6384250122390625, 2.5426674762085515, 2.449262959529336, 2.369463226806091, 2.3122892521514444, 2.283050333978335, 2.270877736449921, 2.2661398792013117, 2.2646365007195577, 2.264415072441482, 2.264172966387209, 2.2641556432663523, 2.2641526812417814, 2.264398548082754, 2.2649516200713413, 2.265520717755241, 2.265920517352683, 2.2661124220075326, 2.2662022870874168, 2.2664171997500953, 2.2666771940569213, 2.266955430531363, 2.267195065288975, 2.2673440197309276, 2.2673518222571007, 2.2673308247750157, 2.2673009131470803, 2.2672482338816566, 2.2670670984298362, 2.2666922156326996, 2.2661636492018538, 2.2657337051235316, 2.265553842004543, 2.2655307159180857, 2.2655938039722843, 2.26579494688254, 2.266395591833907, 2.267530315673026, 2.2687399432098094, 2.269474559281484, 2.2696858933696826, 2.2696847800188578, 2.269487150649506, 2.2686883572651593, 2.2660458086593165, 2.2582631088265734, 2.236502991095803, 2.190955403621243, 2.1241356754780165, 2.043066614858798, 1.9569333851373547, 1.8758643245184075, 1.8090445963753141, 1.763497008900806, 1.7417368911722835, 1.733954191340176, 1.7313116427352884, 1.7305128493515511, 1.730315219981628, 1.7303141066302856, 1.7305254407194954, 1.7312600567917453, 1.7324696843283347, 1.73360440816738, 1.7342050531170665, 1.7344061960270123, 1.7344692840812908, 1.7344461579948052, 1.7342662948774115, 1.7338363508008645, 1.7333077843696465, 1.7329329015718171, 1.7327517661199126, 1.7326990868549703, 1.7326691752271244, 1.732648177744926, 1.7326559802694121, 1.7328049347085135, 1.7330445694657173, 1.733322805942058, 1.733582800250575, 1.7337977129176105, 1.733887577998848, 1.7340794826498784, 1.7344792822432307, 1.735048379925165, 1.735601451915893, 1.7358473187607164, 1.7358443567385913, 1.735827033617042, 1.7355849275537545, 1.7353634992678093, 1.7338601207925701, 1.7291222635625438, 1.716949666016267, 1.6877107478488225, 1.630536773193283, 1.5507370404733194, 1.4573325237925119, 1.3615749877623482, 1.2746852988347757, 1.205713435198065, 1.1611953063895466, 1.1415281967726116, 1.1353201488838707, 1.133024087681783, 1.1319138497096726, 1.1313125630112415, 1.1311042385200525, 1.131217624846748, 1.1315390414667301, 1.1320819577037695, 1.1326782618214672, 1.1332161333080584, 1.1335557439268236, 1.1336470609940874, 1.1336310630879893, 1.133459682047254, 1.1331059356399402, 1.1327429942848692, 1.1325143690162742, 1.1324531083818512, 1.1324429008827455, 1.1325453916767025, 1.1327821224846675, 1.1330011550691792, 1.1331062637904448, 1.1331289140523757, 1.1331262041387824, 1.1330221041773214, 1.1328626089901659, 1.1327416952569926, 1.1327151239571922, 1.1327113900720063, 1.1327468325789778, 1.1328533574336883, 1.1329733954138916, 1.1330279098249783, 1.1330333035797402, 1.1330250688226426, 1.1329640640906482, 1.132870097642622, 1.1327870199911532, 1.1327692808848875, 1.1327707951176493, 1.132797902485734, 1.1328723273220924, 1.1329712523217748, 1.1330418700479132, 1.1330546908913088, 1.1330528304037712, 1.133029715730126, 1.132974944167338, 1.1329077627522597, 1.1328672621431424, 1.1328615879630164, 1.1328636514612385, 1.1328851978086538, 1.1329342646550424, 1.1330016662035525, 1.1330647233263678, 1.1331096156320772, 1.1331433281508334, 1.1331813978585408, 1.133236416218356, 1.1333035882425722, 1.1333562675031317, 1.133370773162378, 1.1333701559644191, 1.1333508513964206, 1.1332432551191398, 1.132898921925896, 1.131966203563651, 1.1289384619946887, 1.1218927980085103, 1.109386903193744, 1.0921174690619964, 1.0722468393446496, 1.051827388407045, 1.0330929938990197, 1.0178268959876138, 1.007361562677124, 1.0023250867817173, 1.000716106832348, 1.0002203396724485, 1.000067697569734, 1.0000207725461103, 1.0000063663960157, 1.0000019486550318, 1.0000005955992908, 1.0000001817562627, 1.0000000553697603, 1.0000000168355723, 1.0000000051082387, 1.0000000015463644, 1.0000000004669227, 1.0000000001406277, 1.0000000000421079, 1.0000000000129612, 1.0000000000028355];
w0x = [6.047425329335867*10^-12, 2.7622598326205093*10^-11, 8.971945566619242*10^-11, 2.996218455513306*10^-10, 9.947911373773268*10^-10, 3.294555768902938*10^-9, 1.088316171380393*10^-8, 3.586835601115286*10^-8, 1.1796582553377614*10^-7, 3.8723349073042977*10^-7, 1.2689301333225753*10^-6, 4.151628651242914*10^-6, 0.000013563668880243794, 0.00004425611320429948, 0.00014423033336540588, 0.00046943583146414757, 0.001525672606169036, 0.004953620155431212, 0.01568388130061356, 0.037980376473787325, 0.07050494757211306, 0.11041875734713621, 0.15392259433760183, 0.1962571604038163, 0.23304992037351555, 0.25969389446723645, 0.2747047212363091, 0.281154934562165, 0.2831420109001687, 0.28387563021860307, 0.28410510965389235, 0.2841466107360923, 0.284148057310035, 0.28411728643307566, 0.28400520771758464, 0.28386215501034073, 0.2837448667345721, 0.2836635815689307, 0.2835914829373022, 0.2834954951521297, 0.28336051021149, 0.2832161482073206, 0.28311081901923907, 0.2830655598375017, 0.28306161830128973, 0.28307450600603057, 0.28316238383604503, 0.28330750400295557, 0.2834263054117746, 0.28347497716486925, 0.28347829902736066, 0.2834498210660966, 0.2832969607643648, 0.28308288267156334, 0.28292013693786083, 0.282861077472262, 0.2828578404768831, 0.2828965873100151, 0.2830762103154318, 0.28328146640242036, 0.28341910500173473, 0.2834430918832124, 0.28343375471064247, 0.2833215432178755, 0.28306480809293494, 0.2828300630367386, 0.2827281555184525, 0.28272607306615927, 0.2827710040004432, 0.2830107943898573, 0.2833418714198501, 0.2835795957485777, 0.28362781575896306, 0.28359850083807725, 0.28344370114268014, 0.28302589646808063, 0.2825874933692611, 0.28234491908942, 0.28232382631576924, 0.28240867528797575, 0.28278017772814434, 0.28338995625810204, 0.28393725785966273, 0.28419767064461204, 0.28421961630004017, 0.2840494227364603, 0.28346743362223425, 0.2826482995526244, 0.28187503886633686, 0.28134683713111874, 0.2811998880610453, 0.2813351320521148, 0.2817935189286057, 0.2829034197732935, 0.2843801927448206, 0.2859149383242225, 0.2898222347853528, 0.305216820972383, 0.34225902546210296, 0.398315603896769, 0.46889751435989124, 0.546801369839831, 0.6227396200256206, 0.6877246441016409, 0.7341141355352055, 0.7590692793593368, 0.767985486412926, 0.7705097919488196, 0.7725768281945382, 0.7738763991907297, 0.7740975057670959, 0.7737420961805671, 0.7731760065840367, 0.7726628489054622, 0.7723057743254721, 0.7721582793689126, 0.7719941278400071, 0.7715780129142584, 0.7712051177745487, 0.7710714784048492, 0.7711183672854439, 0.7712160398176792, 0.7712812635368804, 0.7711420481168525, 0.7709164348207809, 0.7707860973961846, 0.7708227520959186, 0.7709820544634224, 0.771168511787555, 0.7713331669436895, 0.7713582977005319, 0.7713233046852281, 0.7712892282179602, 0.7712876252399424, 0.7712401298771393, 0.7711915630593262, 0.7711907595571144, 0.7712031136136555, 0.7712484530675419, 0.7713029491117014, 0.7713286825515012, 0.7712905246503493, 0.7712611642472622, 0.7711570729692143, 0.7710657481788575, 0.771031233522174, 0.771089202399877, 0.7712437696225873, 0.7713960426224872, 0.7714793549359722, 0.7714793549349688, 0.7713960426201202, 0.7712437696198835, 0.7710892023978188, 0.7710312335239984, 0.7710657481820103, 0.771157072973529, 0.7712611642512227, 0.7712905246519499, 0.7713286825513257, 0.7713029491082118, 0.7712484530604118, 0.7712031136068427, 0.7711907595515002, 0.7711915630566514, 0.7712401298769666, 0.7712876252395814, 0.7712892282177374, 0.7713233046822869, 0.7713582976983491, 0.771333166941611, 0.7711685117870386, 0.7709820544647099, 0.7708227520974802, 0.7707860973968471, 0.7709164348211316, 0.771142048118326, 0.7712812635425174, 0.7712160398241612, 0.7711183672890405, 0.7710714784038476, 0.771205117767404, 0.7715780129054182, 0.7719941278383233, 0.77215827937723, 0.7723057743338401, 0.7726628489045093, 0.7731760065785318, 0.7737420961692503, 0.7740975057599375, 0.7738763991988821, 0.7725768282146669, 0.7705097919516528, 0.7679854863846836, 0.7590692793562366, 0.7341141355276436, 0.6877246440967364, 0.6227396200304566, 0.5468013698439018, 0.4688975143638768, 0.39831560389595183, 0.3422590254556141, 0.30521682096822755, 0.289822234791531, 0.28591493833099, 0.28438019274772486, 0.28290341977237304, 0.2817935189262773, 0.28133513204981375, 0.28119988805907953, 0.2813468371296976, 0.28187503886683973, 0.2826482995546732, 0.28346743362499016, 0.2840494227367152, 0.2842196162996472, 0.28419767064167234, 0.2839372578530093, 0.28338995625125407, 0.2827801777231201, 0.2824086752891091, 0.2823238263194051, 0.28234491909965964, 0.28258749338206707, 0.2830258964797382, 0.2834437011444273, 0.28359850083787724, 0.2836278157536717, 0.2835795957464937, 0.283341871426273, 0.2830107944014901, 0.28277100400612115, 0.2827260730711576, 0.2827281555086843, 0.2828300630223068, 0.2830648080839545, 0.2833215432193714, 0.28343375471333243, 0.28344309189038, 0.28341910500445816, 0.28328146640689855, 0.2830762103205672, 0.2828965873197546, 0.28285784048666285, 0.28286107747847744, 0.28292013693144774, 0.28308288266290904, 0.2832969607552441, 0.2834498210653018, 0.2834782990277116, 0.28347497716749415, 0.28342630540959857, 0.2833075039984061, 0.28316238383406284, 0.2830745060095161, 0.2830616183054649, 0.2830655598428913, 0.28311081901936985, 0.2832161482033649, 0.28336051020389463, 0.2834954951407528, 0.28359148292551484, 0.2836635815627927, 0.2837448667467808, 0.2838621550288979, 0.2840052077316453, 0.28411728641976897, 0.2841480572926249, 0.28414661070681535, 0.2841051096759237, 0.2838756302336005, 0.2831420108961783, 0.28115493454470886, 0.27470472122536027, 0.25969389446784186, 0.23304992038536038, 0.19625716041385477, 0.15392259434230054, 0.11041875733971832, 0.07050494757625064, 0.037980376458370824, 0.015683881323428964, 0.004953620144305032, 0.0015256726104484876, 0.00046943583015493333, 0.00014423033355648708, 0.00004425611333179645, 0.000013563668759207526, 4.151628706303872*10^-6, 1.2689300946688313*10^-6, 3.872335035203604*10^-7, 1.1796582098789762*10^-7, 3.5868352313056926*10^-8, 1.0883157523334361*10^-8, 3.294546813679733*10^-9, 9.947857878197327*10^-10, 2.996107818614762*10^-10, 8.971335691789386*10^-11, 2.7615856405160902*10^-11, 6.042414052720739*10^-12]
w0y =[0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.];
w1 = [4.319032583819702*10^-12, 1.9728793587645348*10^-11, 6.408082671543841*10^-11, 2.1400273790931265*10^-10, 7.105241140154889*10^-10, 2.3531218184480682*10^-9, 7.773244002666611*10^-9, 2.5618810448384393*10^-8, 8.425652668085732*10^-8, 2.765797225590971*10^-7, 9.063273770028556*10^-7, 2.965281349306296*10^-6, 9.68778698151652*10^-6, 0.00003160972131575593, 0.00010301584114592806, 0.00033529234062598353, 0.0010897045215099075, 0.003538099977830172, 0.011202138809018663, 0.027127306103997932, 0.0503578288203912, 0.07886607951731503, 0.10993848856503148, 0.14017574952262296, 0.16645486270562634, 0.18548523576254505, 0.1962066214365551, 0.20081327011626612, 0.20223244879606844, 0.20275644190398226, 0.2029205792507352, 0.20295055324948588, 0.20295171377775334, 0.20292985953347942, 0.20284994970783887, 0.20274782856678453, 0.20266399122131756, 0.20260577268229768, 0.20255402401121214, 0.2024851468075882, 0.20238814187758714, 0.20228432620086306, 0.20220836496685457, 0.2021766459161205, 0.20217424974972653, 0.20218419649493174, 0.20224844858820007, 0.20235389683310095, 0.20244065808706677, 0.2024749991442546, 0.202476849515755, 0.2024554102888114, 0.20234405819440132, 0.202187785192124, 0.20206784073856826, 0.2020243395026846, 0.2020221113653797, 0.20205066260326943, 0.20218142738421335, 0.2023326346038813, 0.20243809682082495, 0.2024612067678236, 0.20245668293155636, 0.20237940217671038, 0.20219483311014688, 0.20202135108134667, 0.20192503928547884, 0.201912339717171, 0.20193411385509547, 0.2020887764576929, 0.20231617594050377, 0.20249600735610443, 0.2025768636785381, 0.2025856070770007, 0.20250878271829795, 0.20227311709020257, 0.2020154698457882, 0.20184061683937843, 0.2017732777503581, 0.20178573028967772, 0.20194866641935078, 0.20222843010337907, 0.20244113304311478, 0.20252685833171225, 0.20252277413118644, 0.20241773922205586, 0.20210817900448513, 0.201870186878932, 0.20180886681428992, 0.20199031153044267, 0.2023102200998205, 0.20270289687767487, 0.20327667890729303, 0.20360178112749505, 0.20348465551149222, 0.20178950475224028, 0.19607498451670535, 0.18273152031662354, 0.15622459302763006, 0.11265398083684719, 0.057442222964292224, -0.0031524368863517734, -0.06226522660943873, -0.11292904004278262, -0.14928523834045002, -0.16841347482971128, -0.17526130793484374, -0.17731893549925448, -0.17882704922168624, -0.17979729205275044, -0.17995990211952412, -0.17975756272827403, -0.1793431544292956, -0.1788786530135476, -0.17858377056337832, -0.17840405479256383, -0.17821827301716664, -0.17786620356879837, -0.1775597890392126, -0.1774367110657056, -0.1774926275938032, -0.17767115587123883, -0.17781532048522625, -0.1778891800461686, -0.177887652997859, -0.17787979038917023, -0.1779772590233493, -0.17805100063155205, -0.1780339593498509, -0.17787864938625164, -0.17751804642562763, -0.17715303287097683, -0.1769653186289765, -0.17692697210993735, -0.1769484495429121, -0.17704867312344302, -0.17742736014477456, -0.17820086655697898, -0.1790475937311529, -0.17962162358815575, -0.1798200928453454, -0.17984187056652182, -0.1797101856381096, -0.17914936966803968, -0.17734023308623434, -0.17209925149453423, -0.1575621430696434, -0.12720435636031663, -0.08269751795516202, -0.028694086034572734, 0.02869408604120119, 0.08269751796130692, 0.12720435636554647, 0.15756214307377314, 0.17209925149608749, 0.17734023308753477, 0.1791493696680376, 0.17971018563703345, 0.17984187056533515, 0.17982009284460715, 0.1796216235854603, 0.1790475937269692, 0.17820086655321465, 0.17742736014166893, 0.17704867312253691, 0.17694844954411137, 0.17692697211180616, 0.17696531863140136, 0.17715303287066148, 0.17751804642307406, 0.17787864938382442, 0.17803395934809962, 0.17805100062997, 0.17797725902219463, 0.17787979038838614, 0.1778876529970052, 0.17788918004663484, 0.1778153204865663, 0.1776711558736929, 0.17749262759365378, 0.17743671106359155, 0.17755978903227942, 0.17786620356140392, 0.1782182730158475, 0.17840405479668, 0.17858377056772584, 0.17887865301279932, 0.17934315442437362, 0.17975756272098492, 0.17995990211580914, 0.17979729206124231, 0.17882704924006948, 0.17731893550535524, 0.17526130791523462, 0.16841347483394714, 0.1492852383343674, 0.11292904004054591, 0.0622652266069956, 0.0031524368887724994, -0.057442222961443024, -0.11265398083441455, -0.1562245930263802, -0.18273152031450773, -0.19607498452329092, -0.20178950475786475, -0.2034846555141398, -0.2036017811274978, -0.20327667890542628, -0.20270289687702556, -0.20231022010050698, -0.20199031153081012, -0.20180886681449192, -0.20187018687820918, -0.20210817900353872, -0.20241773922147646, -0.20252277413045594, -0.2025268583322644, -0.2024411330450718, -0.2022284301064384, -0.2019486664211086, -0.2017857302904637, -0.2017732777483197, -0.2018406168327655, -0.20201546983864963, -0.20227311708168227, -0.2025087827170914, -0.20258560707567438, -0.2025768636829197, -0.20249600735845089, -0.20231617593642, -0.20208877645001624, -0.2019341138506093, -0.2019123397130414, -0.20192503929192954, -0.2020213510916039, -0.20219483311602157, -0.20237940217553174, -0.2024566829297982, -0.20246120676219392, -0.20243809681861455, -0.20233263460080667, -0.20218142738193978, -0.2020506625961071, -0.2020221113577447, -0.2020243394977832, -0.20206784074310716, -0.20218778519771524, -0.20234405820117407, -0.20245541028996436, -0.20247684951619355, -0.20247499914301165, -0.20244065808963874, -0.20235389683714922, -0.20224844858963908, -0.202184196491206, -0.20217424974512538, -0.20217664591052747, -0.2022083649650972, -0.202284326202708, -0.2023881418828237, -0.20248514681634955, -0.20255402402076267, -0.20260577268818009, -0.202663991214249, -0.2027478285549281, -0.2028499496984075, -0.20292985954249987, -0.2029517137890829, -0.20295055326920822, -0.20292057923357354, -0.20275644189153588, -0.2022324487979734, -0.20081327012608802, -0.19620662144187853, -0.18548523575965997, -0.16645486269756185, -0.14017574951918804, -0.10993848856625849, -0.07886607952585803, -0.05035782881852068, -0.02712730611459484, -0.011202138792038074, -0.0035380999855094, -0.0010897045183312362, -0.0003352923415512853, -0.00010301584104271678, -0.000031609721246214474, -9.687787061469601*10^-6, -2.9652812945751857*10^-6, -9.063273957711888*10^-7, -2.765797132455531*1010^-7, -8.425653415073026*10^-8, -2.5618804554513212*10^-8, -7.773243895566488*10^-9, -2.353115364298743*10^-9, -7.105216096100281*10^-10, -2.139959058298512*10^-10, -6.4077757209275*10^-11, -1.9724711694576274*10^-11, -4.315858834561008*10^-12]
w0xx = [-2.863068341270974*10^-12, -1.3077700352939723*10^-11, -4.247706568498021*10^-11, -1.4185471259642427*10^-10, -4.709800841554107*10^-10, -1.5597966072659029*10^-9, -5.152593643602592*10^-9, -1.6981755422008683*10^-8, -5.585051379442783*10^-8, -1.8333439775975804*10^-7, -6.007706633031882*10^-7, -1.965574567967411*10^-6, -6.4216731786863565*10^-6, -0.000020952906989506004, -0.00006828536437940392, -0.00022225280498782004, -0.0007223245423800412, -0.002345274701802667, -0.007425480598920249, -0.01798168083386021, -0.03338032908829769, -0.05227738716759709, -0.07287413116315282, -0.09291728743134231, -0.11033670636791962, -0.12295121355756186, -0.13005803343472128, -0.13311178922703734, -0.13405254782448214, -0.13439987857731553, -0.13450857258944118, -0.13452828433483724, -0.13452899545805946, -0.13451445138021678, -0.13446141536666667, -0.13439369745379004, -0.13433815461306226, -0.1342996391061585, -0.13426545486311314, -0.13421994755436178, -0.13415592377200442, -0.13408743894450942, -0.13403742758236473, -0.13401611364870802, -0.13401432514259137, -0.13402057579028748, -0.13406248731585868, -0.13413155238354946, -0.1341881528222345, -0.13421110861422614, -0.13421258403261685, -0.13419888952424458, -0.13412611575886232, -0.13402403382948586, -0.1339462908713659, -0.1339180782289361, -0.133916594691419, -0.13393510307995937, -0.13402062041384172, -0.13411865886901284, -0.13418524101854867, -0.13419774687830868, -0.13419380099512063, -0.13414130272028196, -0.1340195750370145, -0.13390689482409684, -0.13385431456049465, -0.13385105369990502, -0.13387049725580444, -0.13398111472540206, -0.13413513732259763, -0.13425129445350126, -0.13428174645836644, -0.13427495747038387, -0.13420554498850956, -0.13402171106415012, -0.13382620833208633, -0.13370716589320303, -0.1336887786423064, -0.13372111240136123, -0.13387468587285592, -0.13413443445747203, -0.13435913351708895, -0.1344641515235524, -0.13447085715843723, -0.13439369173114252, -0.13413776564745178, -0.1338241789699861, -0.13352686159737184, -0.13341064907578487, -0.13344614209081973, -0.13352167969659892, -0.13379231878061726, -0.1342643882757578, -0.13483565783325374, -0.1350303676595977, -0.13512663127489455, -0.13763864764777725, -0.14477783743889852, -0.15488178594947705, -0.16759098329905342, -0.18169330536026312, -0.19541027536286787, -0.20717802108304373, -0.21574793760999483, -0.21963110835830005, -0.22178542342505037, -0.22288378775961132, -0.2228432556421555, -0.22241753149024912, -0.22243417620504857, -0.22259500172695576, -0.22269151696789063, -0.22268815265999997, -0.2224984709915089, -0.22221380296636797, -0.22204606834065096, -0.22204330223495988, -0.22212299376827624, -0.22214227157230354, -0.22213118057110193, -0.22212586851427474, -0.2222149836248316, -0.2224084882790261, -0.22260673901903644, -0.2227085578150723, -0.22271898383787042, -0.22264375809054135, -0.2224600551211503, -0.222144832701149, -0.2217585284282061, -0.22143158736811566, -0.22128381396902091, -0.22124872173429033, -0.22125495262698738, -0.22138511602581817, -0.22187490773397944, -0.22285315165380032, -0.223906953381584, -0.2245322570708243, -0.22473517913369856, -0.2247531839617619, -0.2246185878090648, -0.22395096815582893, -0.22174167016812915, -0.2152298590124839, -0.19707676948863845, -0.15908493393989273, -0.10339242970933624, -0.03584955483090774, 0.035849554826342084, 0.10339242970523224, 0.15908493393643008, 0.19707676948569108, 0.21522985901141994, 0.22174167016743104, 0.22395096815608104, 0.2246185878101264, 0.22475318396284022, 0.22473517913436292, 0.22453225707233465, 0.22390695338363129, 0.22285315165574743, 0.22187490773590088, 0.22138511602608452, 0.22125495262613476, 0.22124872173317814, 0.22128381396768798, 0.22143158736851626, 0.2217585284301569, 0.22214483270288737, 0.2224600551225674, 0.22264375809111148, 0.22271898383806535, 0.22270855781541254, 0.22260673901949618, 0.22240848827858772, 0.22221498362303116, 0.22212586851237576, 0.22213118057081016, 0.22214227157352212, 0.22212299377297554, 0.22204330223853375, 0.22204606834063012, 0.22221380296399013, 0.2224984709871271, 0.22268815266077743, 0.22269151697165934, 0.22259500173353453, 0.22243417620896752, 0.2224175314818987, 0.22284325563360535, 0.2228837877634912, 0.2217854234362405, 0.21963110834861876, 0.21574793761467548, 0.20717802108510924, 0.19541027536499567, 0.18169330535964323, 0.16759098329762398, 0.15488178594700602, 0.1447778374363637, 0.13763864764638897, 0.13512663127971986, 0.13503036766418514, 0.13483565783496027, 0.1342643882750016, 0.13379231877955383, 0.1335216796962082, 0.13344614209080596, 0.13341064907605857, 0.13352686159785526, 0.13382417896937887, 0.1341377656463295, 0.13439369173033008, 0.13447085715748908, 0.13446415152379607, 0.13435913351894585, 0.13413443445994788, 0.13387468587518658, 0.13372111240160578, 0.13368877864140105, 0.1337071658893414, 0.133826208326896, 0.13402171105868727, 0.13420554498818543, 0.134274957470574, 0.13428174646167568, 0.13425129445491796, 0.13413513731964186, 0.1339811147200728, 0.13387049725237848, 0.1338510536965422, 0.13385431456424954, 0.13390689483063192, 0.1340195750414211, 0.1341413027200968, 0.13419380099436426, 0.1341977468750249, 0.13418524101763088, 0.13411865886742289, 0.13402062041247995, 0.13393510307514878, 0.13391659468605155, 0.1339180782252387, 0.1339462908736873, 0.1340240338332399, 0.13412611576324976, 0.13419888952476217, 0.1342125840325244, 0.13421110861297303, 0.13418815282360907, 0.1341315523861217, 0.13406248731673742, 0.13402057578780874, 0.13401432513961828, 0.13401611364518615, 0.13403742758146653, 0.1340874389459255, 0.13415592377562374, 0.134219947560277, 0.13426545486950744, 0.13429963911006304, 0.13433815460832937, 0.134393697445862, 0.13446141536035133, 0.1345144513861454, 0.13452899546552652, 0.1345282843478813, 0.1345085725780326, 0.13439987856902644, 0.13405254782573392, 0.13311178923357989, 0.13005803343830782, 0.12295121355566548, 0.11033670636257426, 0.09291728742906061, 0.0728741311639643, 0.052277387173258694, 0.03338032908705715, 0.017981680840884036, 0.007425480587664337, 0.0023452747068916334, 0.0007223245402727493, 0.00022225280560084706, 0.0000682853643113738, 0.000020952906943274732, 6.421673231461814*10^-6, 1.9655745312021536*10^-6, 6.007706753960746*10^-7, 1.8333439108601575*10^-7, 5.585051847510695*10^-8, 1.6981751456627272*10^-8, 5.152593840208675*10^-9, 1.5597924903389668*10^-9, 4.709782418345391*10^-10, 1.4184977563831234*10^-10, 4.247462523084614*10^-11, 1.3074664333251657*10^-11, 2.860754829664413*10^-12]

pd = PlotData1D(sol; solution_variables=cons2cons)

# w0 plot
plot(pd.x, pd.data[:,1], xlims = (-1.0,1.0), label = "DGSEM", title ="w0, t = 0.4", seriestype = :scatter, markersize=1)
plot!(x,w0, label = "FVV")

# w0x plot
plot(pd.x, pd.data[:,2], xlims = (-1.0,1.0), label = "DGSEM", title ="w0x, t = 0.4", seriestype = :scatter, markersize=1)
plot!(x,w0x, label = "FVV")

# w0y plot
plot(pd.x, pd.data[:,3], xlims = (-1.0,1.0), label = "DGSEM", title ="w0y, t = 0.4", seriestype = :scatter, markersize=1)
plot!(x,w0y, label = "FVV")

# w0xy plot 
#plot(pd.x, pd.data[:,7], xlims = (-1.0,1.0), label = "DGSEM", title ="w0y, t = 0.4", seriestype = :scatter, markersize=1)


#savefig("plot0,05.pdf")

# #plot(pd)
# #plot!(getmesh(pd))


