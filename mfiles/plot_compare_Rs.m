%plot
semilogy(saturation_NSL_ay,surface_resistance_ay);
hold on
semilogy(saturation_NSL_ay,surface_resistance_van);
semilogy(saturation_NSL_ay,surface_resistance_Sch);
semilogy(saturation_NSL_cm,surface_resistance_cm);

legend ('new', 'van', 'Sch','cm')

xlim([0 1.01])
ylim([5e-1 10e3])