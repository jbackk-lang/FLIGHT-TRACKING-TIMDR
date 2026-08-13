from timdr_flight import TIMDRFlight

timdr = TIMDRFlight()

# przykładowy tor lotu: [lat, lon, alt, t]
track = [
    [50.0, 19.9, 1000, 0],
    [50.01, 19.91, 1200, 10],
    [50.03, 19.93, 1500, 20],
    [50.06, 19.96, 1800, 30],
    [50.10, 20.00, 2000, 40],
]

flow = timdr.timdr_flow(track)
twist = timdr.twist(track)
stable = timdr.trm_reduce(track)
pred = timdr.predict(track)
diag = timdr.diagnostics(track)

print("TIMDR-flow:\n", flow)
print("Twist:", twist)
print("Stabilized track:\n", stable)
print("Prediction [lat, lon, alt]:\n", pred)
print("Kurs (deg):", diag["course_deg"].round(1))
print("Predkosc wzgledem ziemi (kt):", diag["ground_speed_kt"].round(1))
print("Predkosc pionowa (ft/min):", diag["climb_rate_fpm"].round(0))
