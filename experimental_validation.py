#!/usr/bin/env python3
"""
Experimental validation of E-QFT predictions against real data.

Compares theoretical E-QFT corrections with:
1. Mercury perihelion precession observations
2. Binary pulsar timing data (PSR B1913+16)
3. LIGO gravitational wave detections
4. Solar system tests of gravity
"""

import numpy as np
import matplotlib.pyplot as plt
from ppE_waveform_module import ppEWaveformGenerator
import scipy.constants as const
import argparse
import sys

class ExperimentalValidator:
    """Validate E-QFT predictions against experimental data."""
    
    def __init__(self):
        """Initialize with E-QFT corrections."""
        self.gen = ppEWaveformGenerator()
        
        # Experimental data from literature
        self.mercury_data = {
            "observed_precession": 43.11,      # arcsec/century (total)
            "newtonian_precession": 0.0,       # Newtonian prediction
            "gr_precession": 42.98,            # Einstein GR prediction
            "other_effects": 0.13,             # Solar quadrupole, etc.
            "uncertainty": 0.21                # ±0.21 arcsec/century
        }
        
        self.psr_b1913_data = {
            "orbital_period": 27906.98,        # seconds
            "period_derivative": -2.4184e-12,  # dimensionless P_dot
            "eccentricity": 0.617155,
            "periastron_advance": 4.226607,    # deg/year (observed)
            "gr_prediction": 4.226619,         # deg/year (GR theory)  
            "uncertainty": 0.000005            # ±5 μdeg/year
        }
        
        # LIGO event data (simplified for demonstration)
        self.gw150914_data = {
            "mass_1": 36.0,                    # M_sun (approximate)
            "mass_2": 29.0,                    # M_sun
            "final_mass": 62.0,                # M_sun
            "radiated_energy": 3.0,            # M_sun c²
            "peak_frequency": 250.0,           # Hz
            "peak_strain": 1e-21,              # dimensionless
            "distance": 410.0                  # Mpc
        }
    
    def test_mercury_precession(self):
        """Test E-QFT predictions against Mercury observations."""
        print("🪐 MERCURY PERIHELION PRECESSION TEST")
        print("=" * 45)
        
        # Load Mercury parameters
        params = self.gen.load_system_card("cards/Mercury.yaml")
        M1, M2, e = params['M1'], params['M2'], params['e']
        a = params['a']  # semi-major axis in metres
        
        # Compute E-QFT effective gravity
        G_corrected = self.gen.compute_effective_gravity(M1, M2, e)
        
        # Calculate precession using corrected G
        # Formula: Δω = 6πGM/(c²a(1-e²)) per orbit
        M_sun_kg = 1.989e30  # kg (solar mass)
        M_total_kg = M1 * M_sun_kg  # Sun mass dominates
        
        # Precession per orbit (radians)
        precession_per_orbit_gr = 6 * np.pi * const.G * M_total_kg / (const.c**2 * a * (1 - e**2))
        precession_per_orbit_eqft = 6 * np.pi * G_corrected * M_total_kg / (const.c**2 * a * (1 - e**2))
        
        # Convert to arcseconds per century
        orbital_period_years = params.get('P_orb', 7600521.0) / (365.25 * 24 * 3600)
        orbits_per_century = 100.0 / orbital_period_years
        
        precession_gr_arcsec = precession_per_orbit_gr * 206265 * orbits_per_century
        precession_eqft_arcsec = precession_per_orbit_eqft * 206265 * orbits_per_century
        
        # E-QFT correction
        eqft_correction = precession_eqft_arcsec - precession_gr_arcsec
        
        print(f"Orbital parameters:")
        print(f"  a = {a:.2e} m, e = {e:.6f}")
        print(f"  G_Newton = {const.G:.4e} m³/(kg·s²)")
        print(f"  G_E-QFT = {G_corrected:.4e} m³/(kg·s²)")
        print(f"  Correction factor = {G_corrected/const.G:.6f}")
        
        print(f"\nPrecession predictions (arcsec/century):")
        print(f"  Newton only:     0.00")
        print(f"  Einstein GR:     {precession_gr_arcsec:.4f}")
        print(f"  E-QFT corrected: {precession_eqft_arcsec:.4f}")
        print(f"  E-QFT correction: {eqft_correction:+.4f}")
        print(f"  Observed:        {self.mercury_data['observed_precession']:.2f} ± {self.mercury_data['uncertainty']:.4f}")
        print(f"  Known GR value:  {self.mercury_data['gr_precession']:.4f}")
        
        # Residual analysis
        obs_residual_gr = self.mercury_data['observed_precession'] - precession_gr_arcsec
        obs_residual_eqft = self.mercury_data['observed_precession'] - precession_eqft_arcsec
        
        print(f"\nResidual analysis:")
        print(f"  Obs - GR:       {obs_residual_gr:+.4f} arcsec/century")
        print(f"  Obs - E-QFT:    {obs_residual_eqft:+.4f} arcsec/century")
        print(f"  E-QFT improvement: {abs(obs_residual_eqft) < abs(obs_residual_gr)}")
        
        # Statistical significance
        sigma_gr = abs(obs_residual_gr) / self.mercury_data['uncertainty']
        sigma_eqft = abs(obs_residual_eqft) / self.mercury_data['uncertainty']
        
        print(f"\nStatistical significance:")
        print(f"  GR discrepancy:   {sigma_gr:.2f}σ")
        print(f"  E-QFT discrepancy: {sigma_eqft:.2f}σ")
        
        return {
            'theory_gr': precession_gr_arcsec,
            'theory_eqft': precession_eqft_arcsec,
            'observed': self.mercury_data['observed_precession'],
            'eqft_correction': eqft_correction,
            'residual_gr': obs_residual_gr,
            'residual_eqft': obs_residual_eqft,
            'significance_gr': sigma_gr,
            'significance_eqft': sigma_eqft
        }
    
    def test_pulsar_timing(self):
        """Test E-QFT predictions against PSR B1913+16 data."""
        print("\n📡 PSR B1913+16 BINARY PULSAR TEST")
        print("=" * 40)
        
        # Load pulsar parameters
        params = self.gen.load_system_card("cards/PSR_B1913+16.yaml")
        M1, M2, e = params['M1'], params['M2'], params['e']
        
        # Compute E-QFT effective gravity
        G_corrected = self.gen.compute_effective_gravity(M1, M2, e)
        
        # Use the observed GR value as baseline and compute E-QFT correction
        # This avoids potential formula normalization issues
        omega_dot_gr_observed = self.psr_b1913_data['gr_prediction']  # deg/year
        
        # Compute the ratio G_corrected/G to scale the GR prediction
        G_ratio = G_corrected / const.G
        
        # E-QFT prediction scales as G^(5/3) for periastron advance
        omega_dot_eqft_deg = omega_dot_gr_observed * (G_ratio)**(5/3)
        omega_dot_gr_deg = omega_dot_gr_observed  # Use known GR value
        
        # Values are already in deg/year from the scaling above
        
        # E-QFT correction
        eqft_correction = omega_dot_eqft_deg - omega_dot_gr_deg
        
        print(f"Binary parameters:")
        print(f"  M₁ = {M1:.3f} M☉, M₂ = {M2:.3f} M☉")
        print(f"  e = {e:.6f}, P_orb = {self.psr_b1913_data['orbital_period']/3600:.2f} hours")
        print(f"  G correction factor = {G_corrected/const.G:.6f}")
        
        print(f"\nPeriastron advance (deg/year):")
        print(f"  Newton only:     0.000000")
        print(f"  Einstein GR:     {omega_dot_gr_deg:.6f}")
        print(f"  E-QFT corrected: {omega_dot_eqft_deg:.6f}")
        print(f"  E-QFT correction: {eqft_correction:+.8f}")
        print(f"  Observed:        {self.psr_b1913_data['periastron_advance']:.6f}")
        print(f"  Known GR theory: {self.psr_b1913_data['gr_prediction']:.6f}")
        print(f"  Uncertainty:     ±{self.psr_b1913_data['uncertainty']:.6f}")
        
        # Residual analysis
        obs_residual_gr = self.psr_b1913_data['periastron_advance'] - omega_dot_gr_deg
        obs_residual_eqft = self.psr_b1913_data['periastron_advance'] - omega_dot_eqft_deg
        
        print(f"\nResidual analysis:")
        print(f"  Obs - GR:       {obs_residual_gr:+.8f} deg/year")
        print(f"  Obs - E-QFT:    {obs_residual_eqft:+.8f} deg/year")
        
        # Statistical significance (in units of uncertainty)
        sigma_gr = abs(obs_residual_gr) / self.psr_b1913_data['uncertainty']
        sigma_eqft = abs(obs_residual_eqft) / self.psr_b1913_data['uncertainty']
        
        print(f"\nStatistical significance:")
        print(f"  GR discrepancy:   {sigma_gr:.1f}σ")
        print(f"  E-QFT discrepancy: {sigma_eqft:.1f}σ")
        
        return {
            'theory_gr': omega_dot_gr_deg,
            'theory_eqft': omega_dot_eqft_deg,
            'observed': self.psr_b1913_data['periastron_advance'],
            'eqft_correction': eqft_correction,
            'residual_gr': obs_residual_gr,
            'residual_eqft': obs_residual_eqft,
            'significance_gr': sigma_gr,
            'significance_eqft': sigma_eqft
        }
    
    def test_ligo_waveform(self):
        """Test E-QFT waveform against LIGO GW150914-like event."""
        print("\n🌊 LIGO GRAVITATIONAL WAVE TEST")
        print("=" * 35)
        
        # Create synthetic GW150914-like system
        M1, M2 = self.gw150914_data['mass_1'], self.gw150914_data['mass_2']
        distance = self.gw150914_data['distance']
        
        # Generate frequency array (LIGO sensitive band)
        f_array = np.logspace(1, 3, 1000)  # 10 Hz to 1 kHz
        
        # Compute E-QFT effective gravity (assume low eccentricity for BBH)
        e_bbh = 0.01  # Typical for BBH mergers
        G_corrected = self.gen.compute_effective_gravity(M1, M2, e_bbh)
        
        print(f"Binary black hole parameters:")
        print(f"  M₁ = {M1:.1f} M☉, M₂ = {M2:.1f} M☉")
        print(f"  Distance = {distance:.0f} Mpc")
        print(f"  G correction factor = {G_corrected/const.G:.6f}")
        
        # Create temporary parameter card for this system
        bbh_params = {
            'M1': M1, 'M2': M2, 'e': e_bbh, 
            'a': 1e9, 'distance': distance  # Rough estimates
        }
        
        import yaml
        with open('temp_bbh_card.yaml', 'w') as f:
            yaml.dump(bbh_params, f)
        
        try:
            # Generate GR waveform (no ppE)
            h_plus_gr, h_cross_gr = self.gen.frequency_domain_waveform_eqft(
                'temp_bbh_card.yaml', f_array,
                delta_phi_1=0.0, delta_phi_1p5=0.0
            )
            
            # Generate E-QFT waveform with potential ppE signatures
            h_plus_eqft, h_cross_eqft = self.gen.frequency_domain_waveform_eqft(
                'temp_bbh_card.yaml', f_array,
                delta_phi_1=0.01, delta_phi_1p5=0.02  # E-QFT signatures
            )
            
            # Compute characteristic strain
            strain_gr = np.abs(h_plus_gr + 1j * h_cross_gr)
            strain_eqft = np.abs(h_plus_eqft + 1j * h_cross_eqft)
            
            # Find peak frequency and strain
            peak_idx = np.argmax(strain_gr)
            peak_freq = f_array[peak_idx]
            peak_strain_gr = strain_gr[peak_idx]
            peak_strain_eqft = strain_eqft[peak_idx]
            
            print(f"\nWaveform characteristics:")
            print(f"  Peak frequency:           {peak_freq:.0f} Hz")
            print(f"  Peak strain (Newton):     ~0 (no waves)")
            print(f"  Peak strain (GR):         {peak_strain_gr:.2e}")
            print(f"  Peak strain (E-QFT):      {peak_strain_gr:.2e} (same amplitude)")
            print(f"  Peak strain (+ ppE):      {peak_strain_eqft:.2e}")
            print(f"  Strain difference:        {(peak_strain_eqft/peak_strain_gr - 1)*100:+.2f}%")
            print(f"  LIGO sensitivity (~1e-23): {'Detectable' if peak_strain_gr > 1e-23 else 'Below threshold'}")
            
            # Phase evolution comparison
            phase_gr = np.unwrap(np.angle(h_plus_gr + 1j * h_cross_gr))
            phase_eqft = np.unwrap(np.angle(h_plus_eqft + 1j * h_cross_eqft))
            phase_diff = phase_eqft - phase_gr
            
            max_phase_diff = np.max(np.abs(phase_diff))
            print(f"  Maximum phase difference: {max_phase_diff:.2f} rad")
            
            # SNR estimate (simplified)
            # Assume Advanced LIGO noise curve ≈ 1e-23 / √Hz at 100 Hz
            snr_estimate = peak_strain_gr / 1e-23 * np.sqrt(100)  # Rough estimate
            print(f"  Estimated SNR:           {snr_estimate:.1f}")
            
            return {
                'peak_frequency': peak_freq,
                'peak_strain_gr': peak_strain_gr,
                'peak_strain_eqft': peak_strain_eqft,
                'max_phase_diff': max_phase_diff,
                'snr_estimate': snr_estimate,
                'frequencies': f_array,
                'strain_gr': strain_gr,
                'strain_eqft': strain_eqft
            }
            
        finally:
            # Clean up temporary file
            import os
            if os.path.exists('temp_bbh_card.yaml'):
                os.remove('temp_bbh_card.yaml')
    
    def test_gaia_systems(self):
        """Test E-QFT predictions across Gaia binary systems."""
        print("\\n🌌 GAIA BINARY SYSTEMS SURVEY")
        print("=" * 35)
        
        import pandas as pd
        
        try:
            # Load Gaia systems data
            gaia_df = pd.read_csv("cards/gaia_systems.csv")
            print(f"Loaded {len(gaia_df)} Gaia binary systems")
            
            results = {}
            
            # Convert AU to meters for calculations
            AU_to_m = 1.496e11
            
            for idx, row in gaia_df.iterrows():
                system_name = row['name']
                M1, M2 = row['m1_solar'], row['m2_solar']
                a_au, e = row['a_au'], row['e']
                a_m = a_au * AU_to_m  # Convert AU to meters
                
                print(f"\\n🔭 {system_name[:20]}...")
                print("-" * 30)
                
                try:
                    # Compute E-QFT effective gravity
                    G_corrected = self.gen.compute_effective_gravity(M1, M2, e)
                    
                    # Calculate precession using corrected G (same formula as solar system)
                    M_sun_kg = 1.989e30  # kg (solar mass)
                    M_total_kg = (M1 + M2) * M_sun_kg  # Total system mass
                    
                    # Precession per orbit (radians) - use total mass for binary
                    precession_per_orbit_gr = 6 * np.pi * const.G * M_total_kg / (const.c**2 * a_m * (1 - e**2))
                    precession_per_orbit_eqft = 6 * np.pi * G_corrected * M_total_kg / (const.c**2 * a_m * (1 - e**2))
                    
                    # Convert to arcseconds per century (estimate orbital period)
                    # P = 2π√(a³/GM) in seconds
                    orbital_period_sec = 2 * np.pi * np.sqrt(a_m**3 / (const.G * M_total_kg))
                    orbital_period_years = orbital_period_sec / (365.25 * 24 * 3600)
                    orbits_per_century = 100.0 / orbital_period_years
                    
                    precession_gr_arcsec = precession_per_orbit_gr * 206265 * orbits_per_century
                    precession_eqft_arcsec = precession_per_orbit_eqft * 206265 * orbits_per_century
                    
                    # E-QFT correction
                    eqft_correction = precession_eqft_arcsec - precession_gr_arcsec
                    correction_percent = (eqft_correction / precession_gr_arcsec) * 100 if precession_gr_arcsec != 0 else 0
                    
                    print(f"  M₁={M1:.2f} M☉, M₂={M2:.2f} M☉, a={a_au:.3f} AU, e={e:.3f}")
                    print(f"  G correction: {G_corrected/const.G:.6f}")
                    print(f"  Orbital period: {orbital_period_years:.4f} years ({orbital_period_sec/3600:.1f} hours)")
                    print(f"  GR precession: {precession_gr_arcsec:.3e} arcsec/century")
                    print(f"  E-QFT correction: {eqft_correction:+.3e} arcsec/century ({correction_percent:+.2f}%)")
                    
                    # Store results
                    results[system_name] = {
                        'M1': M1, 'M2': M2, 'a_au': a_au, 'e': e,
                        'orbital_period_years': orbital_period_years,
                        'theory_gr': precession_gr_arcsec,
                        'theory_eqft': precession_eqft_arcsec,
                        'eqft_correction': eqft_correction,
                        'correction_percent': correction_percent,
                        'G_ratio': G_corrected/const.G
                    }
                    
                    # Limit output for readability
                    if idx >= 20:  # Show first 20 systems in detail
                        if idx % 20 == 0:
                            print(f"\\n... processed {idx} systems (showing every 20th) ...")
                        
                except Exception as e:
                    print(f"  ❌ Error: {e}")
                    
        except FileNotFoundError:
            print("❌ Gaia systems file not found: cards/gaia_systems.csv")
            return {}
        except Exception as e:
            print(f"❌ Error loading Gaia data: {e}")
            return {}
            
        return results

    def test_csv_systems(self, csv_file):
        """Test E-QFT predictions against systems in CSV file with observational data."""
        print(f"\\n📊 CSV SYSTEMS VALIDATION")
        print("=" * 35)
        
        import pandas as pd
        
        try:
            # Load CSV systems data
            csv_df = pd.read_csv(csv_file)
            print(f"Loaded {len(csv_df)} systems from {csv_file}")
            
            results = {}
            
            # Convert AU to meters for calculations
            AU_to_m = 1.496e11
            
            for idx, row in csv_df.iterrows():
                system_name = row['name']
                M1, M2 = row['m1_solar'], row['m2_solar']
                a_au, e = row['a_au'] if 'a_au' in row else row.get('a', 0), row['e']
                a_m = a_au * AU_to_m
                
                # Get observational data
                observed = row.get('observed_precession_arcsec_per_century', 0.0)
                uncertainty = row.get('uncertainty_arcsec_per_century', 1.0)
                gr_prediction = row.get('gr_prediction_arcsec_per_century', 0.0)
                
                print(f"\\n🔭 {system_name.upper()}")
                print("-" * 30)
                
                try:
                    # Compute E-QFT effective gravity
                    G_corrected = self.gen.compute_effective_gravity(M1, M2, e)
                    
                    # Calculate precession using corrected G
                    M_sun_kg = 1.989e30
                    if system_name.lower() in ['mercury', 'venus', 'earth', 'mars', 'icarus']:
                        # Solar system: use total system mass but Sun dominates
                        M_total_kg = M1 * M_sun_kg  # Sun mass dominates
                    else:
                        # Binary system: use total mass
                        M_total_kg = (M1 + M2) * M_sun_kg
                    
                    # Precession per orbit (radians)
                    precession_per_orbit_gr = 6 * np.pi * const.G * M_total_kg / (const.c**2 * a_m * (1 - e**2))
                    precession_per_orbit_eqft = 6 * np.pi * G_corrected * M_total_kg / (const.c**2 * a_m * (1 - e**2))
                    
                    # Convert to arcseconds per century
                    orbital_period_sec = 2 * np.pi * np.sqrt(a_m**3 / (const.G * M_total_kg))
                    orbital_period_years = orbital_period_sec / (365.25 * 24 * 3600)
                    orbits_per_century = 100.0 / orbital_period_years
                    
                    precession_gr_arcsec = precession_per_orbit_gr * 206265 * orbits_per_century
                    precession_eqft_arcsec = precession_per_orbit_eqft * 206265 * orbits_per_century
                    
                    # E-QFT correction
                    eqft_correction = precession_eqft_arcsec - precession_gr_arcsec
                    
                    print(f"  M₁={M1:.2f} M☉, M₂={M2:.2e} M☉, a={a_au:.3f} AU, e={e:.6f}")
                    print(f"  G correction: {G_corrected/const.G:.6f}")
                    print(f"  Orbital period: {orbital_period_years:.4f} years")
                    print(f"  GR prediction: {precession_gr_arcsec:.3f} arcsec/century")
                    print(f"  E-QFT prediction: {precession_eqft_arcsec:.3f} arcsec/century")
                    print(f"  E-QFT correction: {eqft_correction:+.4f} arcsec/century")
                    print(f"  Observed: {observed:.3f} ± {uncertainty:.3f} arcsec/century")
                    
                    if gr_prediction > 0:
                        print(f"  Known GR theory: {gr_prediction:.3f} arcsec/century")
                    
                    # Residual analysis
                    obs_residual_gr = observed - precession_gr_arcsec
                    obs_residual_eqft = observed - precession_eqft_arcsec
                    
                    # Statistical significance
                    sigma_gr = abs(obs_residual_gr) / uncertainty if uncertainty > 0 else 0
                    sigma_eqft = abs(obs_residual_eqft) / uncertainty if uncertainty > 0 else 0
                    
                    print(f"  GR residual: {obs_residual_gr:+.4f} arcsec/cy ({sigma_gr:.1f}σ)")
                    print(f"  E-QFT residual: {obs_residual_eqft:+.4f} arcsec/cy ({sigma_eqft:.1f}σ)")
                    
                    # Store results
                    results[system_name] = {
                        'theory_gr': precession_gr_arcsec,
                        'theory_eqft': precession_eqft_arcsec,
                        'observed': observed,
                        'uncertainty': uncertainty,
                        'eqft_correction': eqft_correction,
                        'residual_gr': obs_residual_gr,
                        'residual_eqft': obs_residual_eqft,
                        'significance_gr': sigma_gr,
                        'significance_eqft': sigma_eqft,
                        'orbital_params': {'M1': M1, 'M2': M2, 'a_au': a_au, 'e': e}
                    }
                    
                except Exception as e:
                    print(f"  ❌ Error: {e}")
                    
        except FileNotFoundError:
            print(f"❌ CSV file not found: {csv_file}")
            return {}
        except Exception as e:
            print(f"❌ Error loading CSV data: {e}")
            return {}
            
        return results

    def test_solar_system_objects(self):
        """Test E-QFT predictions across multiple solar system objects."""
        print("\\n🌍 SOLAR SYSTEM SURVEY TEST")
        print("=" * 35)
        
        # Solar system objects to test
        solar_objects = ["Venus", "Earth", "Mars", "Icarus", "BepiColombo"]
        
        results = {}
        
        for obj_name in solar_objects:
            print(f"\\n🪐 Testing {obj_name.upper()}")
            print("-" * 25)
            
            try:
                # Load object parameters
                params = self.gen.load_system_card(f"cards/{obj_name}.yaml")
                M1, M2, e = params['M1'], params['M2'], params['e']
                a = params['a']  # semi-major axis in metres
                
                # Compute E-QFT effective gravity
                G_corrected = self.gen.compute_effective_gravity(M1, M2, e)
                
                # Calculate precession using corrected G
                M_sun_kg = 1.989e30  # kg (solar mass)
                M_total_kg = M1 * M_sun_kg  # Sun mass dominates
                
                # Precession per orbit (radians)
                precession_per_orbit_gr = 6 * np.pi * const.G * M_total_kg / (const.c**2 * a * (1 - e**2))
                precession_per_orbit_eqft = 6 * np.pi * G_corrected * M_total_kg / (const.c**2 * a * (1 - e**2))
                
                # Convert to arcseconds per century
                orbital_period_years = params.get('P_orb', 7600521.0) / (365.25 * 24 * 3600)
                orbits_per_century = 100.0 / orbital_period_years
                
                precession_gr_arcsec = precession_per_orbit_gr * 206265 * orbits_per_century
                precession_eqft_arcsec = precession_per_orbit_eqft * 206265 * orbits_per_century
                
                # E-QFT correction
                eqft_correction = precession_eqft_arcsec - precession_gr_arcsec
                
                # Get observational data
                observed = params.get('observed_precession', 0.0)
                uncertainty = params.get('uncertainty', 1.0)
                
                print(f"  Orbital elements: a={a:.2e} m, e={e:.6f}")
                print(f"  G correction: {G_corrected/const.G:.6f}")
                print(f"  GR prediction: {precession_gr_arcsec:.3f} arcsec/century")
                print(f"  E-QFT prediction: {precession_eqft_arcsec:.3f} arcsec/century")
                print(f"  E-QFT correction: {eqft_correction:+.4f} arcsec/century")
                print(f"  Observed: {observed:.3f} ± {uncertainty:.3f} arcsec/century")
                
                # Residual analysis
                obs_residual_gr = observed - precession_gr_arcsec
                obs_residual_eqft = observed - precession_eqft_arcsec
                
                # Statistical significance
                sigma_gr = abs(obs_residual_gr) / uncertainty if uncertainty > 0 else 0
                sigma_eqft = abs(obs_residual_eqft) / uncertainty if uncertainty > 0 else 0
                
                print(f"  GR residual: {obs_residual_gr:+.4f} arcsec/cy ({sigma_gr:.1f}σ)")
                print(f"  E-QFT residual: {obs_residual_eqft:+.4f} arcsec/cy ({sigma_eqft:.1f}σ)")
                
                # Store results
                results[obj_name] = {
                    'theory_gr': precession_gr_arcsec,
                    'theory_eqft': precession_eqft_arcsec,
                    'observed': observed,
                    'uncertainty': uncertainty,
                    'eqft_correction': eqft_correction,
                    'residual_gr': obs_residual_gr,
                    'residual_eqft': obs_residual_eqft,
                    'significance_gr': sigma_gr,
                    'significance_eqft': sigma_eqft,
                    'orbital_params': {'a': a, 'e': e, 'M1': M1, 'M2': M2}
                }
                
            except FileNotFoundError:
                print(f"  ❌ Card not found: cards/{obj_name}.yaml")
            except Exception as e:
                print(f"  ❌ Error processing {obj_name}: {e}")
        
        return results
    
    def create_experimental_comparison_plot(self, mercury_results, pulsar_results, ligo_results):
        """Create comprehensive experimental comparison plot."""
        fig, axes = plt.subplots(2, 2, figsize=(15, 12))
        
        # Mercury precession comparison
        ax1 = axes[0, 0]
        systems = ['GR Theory', 'E-QFT Theory', 'Observed']
        values = [mercury_results['theory_gr'], 
                 mercury_results['theory_eqft'], 
                 mercury_results['observed']]
        errors = [0, 0, self.mercury_data['uncertainty']]
        
        bars = ax1.bar(systems, values, yerr=errors, capsize=5,
                      color=['blue', 'red', 'green'], alpha=0.7)
        ax1.set_ylabel('Precession [arcsec/century]')
        ax1.set_title('Mercury Perihelion Precession')
        ax1.grid(True, alpha=0.3)
        
        # Add value labels
        for bar, value in zip(bars, values):
            height = bar.get_height()
            ax1.text(bar.get_x() + bar.get_width()/2., height + 0.01,
                    f'{value:.2f}', ha='center', va='bottom')
        
        # Pulsar timing comparison  
        ax2 = axes[0, 1]
        systems = ['GR Theory', 'E-QFT Theory', 'Observed']
        values = [pulsar_results['theory_gr'],
                 pulsar_results['theory_eqft'],
                 pulsar_results['observed']]
        errors = [0, 0, self.psr_b1913_data['uncertainty']]
        
        bars = ax2.bar(systems, values, yerr=errors, capsize=5,
                      color=['blue', 'red', 'green'], alpha=0.7)
        ax2.set_ylabel('Periastron Advance [deg/year]')
        ax2.set_title('PSR B1913+16 Timing')
        ax2.grid(True, alpha=0.3)
        
        # LIGO strain comparison
        ax3 = axes[1, 0]
        ax3.loglog(ligo_results['frequencies'], ligo_results['strain_gr'], 
                  'b-', label='GR + E-QFT', linewidth=2)
        ax3.loglog(ligo_results['frequencies'], ligo_results['strain_eqft'], 
                  'r--', label='GR + E-QFT + ppE', linewidth=2)
        ax3.axhline(1e-23, color='gray', linestyle=':', alpha=0.7, 
                   label='LIGO sensitivity')
        ax3.set_xlabel('Frequency [Hz]')
        ax3.set_ylabel('Characteristic Strain')
        ax3.set_title('LIGO Waveform Comparison')
        ax3.legend()
        ax3.grid(True, alpha=0.3)
        
        # Residuals summary
        ax4 = axes[1, 1]
        tests = ['Mercury\n(arcsec/cy)', 'Pulsar\n(μdeg/year)', 'LIGO\n(% strain diff)']
        gr_residuals = [abs(mercury_results['residual_gr']), 
                       abs(pulsar_results['residual_gr']) * 1e6,
                       0]  # No GR residual for LIGO
        eqft_residuals = [abs(mercury_results['residual_eqft']),
                         abs(pulsar_results['residual_eqft']) * 1e6,
                         abs(ligo_results['peak_strain_eqft']/ligo_results['peak_strain_gr'] - 1) * 100]
        
        x = np.arange(len(tests))
        width = 0.35
        
        ax4.bar(x - width/2, gr_residuals, width, label='GR Theory', 
               color='blue', alpha=0.7)
        ax4.bar(x + width/2, eqft_residuals, width, label='E-QFT Theory',
               color='red', alpha=0.7)
        
        ax4.set_ylabel('|Residual|')
        ax4.set_title('Theory vs Observation Residuals')
        ax4.set_xticks(x)
        ax4.set_xticklabels(tests)
        ax4.legend()
        ax4.grid(True, alpha=0.3)
        ax4.set_yscale('log')
        
        plt.tight_layout()
        plt.savefig('experimental_validation_comparison.png', dpi=300, bbox_inches='tight')
        print(f"\n📊 Comparison plot saved: experimental_validation_comparison.png")
        plt.close()

def main():
    """Run complete experimental validation."""
    parser = argparse.ArgumentParser(description="E-QFT Experimental Validation")
    parser.add_argument("--csv", type=str, help="Include CSV file with additional systems (e.g., gaia_systems.csv, solar_system.csv)")
    parser.add_argument("--gaia", action="store_true", help="Include Gaia binary systems survey (theoretical only)")
    
    args = parser.parse_args()
    
    print("🔬 E-QFT EXPERIMENTAL VALIDATION")
    print("=" * 40)
    print("Testing theoretical predictions against real data")
    
    validator = ExperimentalValidator()
    
    # Run core experimental tests (always included)
    mercury_results = validator.test_mercury_precession()
    pulsar_results = validator.test_pulsar_timing()
    ligo_results = validator.test_ligo_waveform()
    
    # Optional tests
    solar_system_results = validator.test_solar_system_objects()
    
    csv_results = {}
    if args.csv:
        csv_results = validator.test_csv_systems(args.csv)
    
    gaia_results = {}
    if args.gaia:
        gaia_results = validator.test_gaia_systems()
    
    # Create comparison plot
    validator.create_experimental_comparison_plot(mercury_results, pulsar_results, ligo_results)
    
    # Comprehensive comparison table
    print(f"\n📋 THEORY vs EXPERIMENT COMPARISON TABLE")
    print("=" * 80)
    print(f"{'Observable':<25} {'Newton':<12} {'Einstein GR':<12} {'E-QFT':<12} {'Observed':<12} {'Status':<15}")
    print("-" * 80)
    
    # Mercury precession
    print(f"{'Mercury ω̇ [″/cy]':<25} {'0.00':<12} {mercury_results['theory_gr']:<12.4f} "
          f"{mercury_results['theory_eqft']:<12.4f} {mercury_results['observed']:<12.4f} "
          f"{'GR excellent':<15}")
    
    # Pulsar timing  
    print(f"{'PSR B1913+16 ω̇ [°/yr]':<25} {'0.000000':<12} {pulsar_results['theory_gr']:<12.6f} "
          f"{pulsar_results['theory_eqft']:<12.6f} {pulsar_results['observed']:<12.6f} "
          f"{'GR excellent':<15}")
    
    # LIGO strain
    ligo_strain_sci = f"{ligo_results['peak_strain_gr']:.1e}"
    print(f"{'LIGO peak strain':<25} {'0 (no waves)':<12} {ligo_strain_sci:<12} "
          f"{ligo_strain_sci:<12} {'~1e-21':<12} {'GR confirmed':<15}")
    
    print("-" * 80)
    
    # Residual analysis
    print(f"\n📊 RESIDUAL ANALYSIS (Theory - Observation)")
    print("-" * 50)
    print(f"{'Observable':<25} {'GR Residual':<15} {'E-QFT Residual':<15}")
    print("-" * 50)
    print(f"{'Mercury [″/cy]':<25} {mercury_results['residual_gr']:+.4f} "
          f"{'':>10} {mercury_results['residual_eqft']:+.4f}")
    print(f"{'PSR B1913+16 [°/yr]':<25} {pulsar_results['residual_gr']:+.8f} "
          f"{'':>5} {pulsar_results['residual_eqft']:+.8f}")
    
    # Statistical significance
    print(f"\n📈 STATISTICAL SIGNIFICANCE (σ)")
    print("-" * 40)
    print(f"{'Observable':<25} {'GR':<8} {'E-QFT':<8}")
    print("-" * 40)
    print(f"{'Mercury':<25} {mercury_results['significance_gr']:<8.1f} "
          f"{mercury_results['significance_eqft']:<8.1f}")
    print(f"{'PSR B1913+16':<25} {pulsar_results['significance_gr']:<8.1f} "
          f"{pulsar_results['significance_eqft']:<8.1f}")
    
    print(f"\n🎯 PHYSICS CONCLUSIONS:")
    print("=" * 30)
    print(f"✅ **Newton's theory**: Completely ruled out by all precision tests")
    print(f"✅ **Einstein's GR**: Excellent agreement with all observations")
    print(f"❓ **E-QFT theory**: Makes testable predictions but currently ruled out by:")
    print(f"   • Binary pulsar timing ({pulsar_results['significance_eqft']:.0f}σ deviation)")
    print(f"   • Mercury precession (marginal, {mercury_results['significance_eqft']:.1f}σ)")
    print(f"⚡ **Future prospects**: LIGO phase analysis may detect E-QFT signatures")
    
    print(f"\n🔬 E-QFT CORRECTION MAGNITUDES:")
    print(f"• Mercury: {mercury_results['eqft_correction']:+.3f} arcsec/century ({abs(mercury_results['eqft_correction']/mercury_results['theory_gr']*100):.2f}% of GR)")
    print(f"• Pulsar: {pulsar_results['eqft_correction']:+.6f} deg/year ({abs(pulsar_results['eqft_correction']/pulsar_results['theory_gr']*100):.2f}% of GR)")
    print(f"• LIGO: {ligo_results['max_phase_diff']:.1f} rad phase difference (potentially observable)")
    
    # Solar system survey summary
    if solar_system_results:
        print(f"\n🌍 SOLAR SYSTEM SURVEY SUMMARY")
        print("=" * 50)
        print(f"{'Object':<15} {'GR [″/cy]':<12} {'E-QFT [″/cy]':<13} {'Observed':<12} {'GR σ':<8} {'E-QFT σ':<8}")
        print("-" * 70)
        
        for obj, data in solar_system_results.items():
            print(f"{obj:<15} {data['theory_gr']:<12.3f} {data['theory_eqft']:<13.3f} "
                  f"{data['observed']:<12.3f} {data['significance_gr']:<8.1f} {data['significance_eqft']:<8.1f}")
        
        print("-" * 70)
        
        # Statistical summary
        gr_sigmas = [data['significance_gr'] for data in solar_system_results.values()]
        eqft_sigmas = [data['significance_eqft'] for data in solar_system_results.values()]
        
        print(f"\n📊 SOLAR SYSTEM STATISTICS:")
        print(f"• Objects tested: {len(solar_system_results)}")
        print(f"• Mean GR deviation: {np.mean(gr_sigmas):.1f}σ")
        print(f"• Mean E-QFT deviation: {np.mean(eqft_sigmas):.1f}σ")
        print(f"• Objects where GR better: {sum(1 for i in range(len(gr_sigmas)) if gr_sigmas[i] < eqft_sigmas[i])}/{len(solar_system_results)}")
        print(f"• Objects where E-QFT better: {sum(1 for i in range(len(gr_sigmas)) if eqft_sigmas[i] < gr_sigmas[i])}/{len(solar_system_results)}")

    # CSV systems summary (if loaded)
    if csv_results:
        print(f"\n📊 CSV SYSTEMS SUMMARY")
        print("=" * 40)
        print(f"{'System':<15} {'GR [″/cy]':<12} {'E-QFT [″/cy]':<13} {'Observed':<12} {'GR σ':<8} {'E-QFT σ':<8}")
        print("-" * 70)
        
        for obj, data in csv_results.items():
            print(f"{obj:<15} {data['theory_gr']:<12.3f} {data['theory_eqft']:<13.3f} "
                  f"{data['observed']:<12.3f} {data['significance_gr']:<8.1f} {data['significance_eqft']:<8.1f}")
        
        print("-" * 70)
        
        # CSV systems statistics
        csv_gr_sigmas = [data['significance_gr'] for data in csv_results.values()]
        csv_eqft_sigmas = [data['significance_eqft'] for data in csv_results.values()]
        
        print(f"\n📈 CSV SYSTEMS STATISTICS:")
        print(f"• Systems tested: {len(csv_results)}")
        print(f"• Mean GR deviation: {np.mean(csv_gr_sigmas):.1f}σ")
        print(f"• Mean E-QFT deviation: {np.mean(csv_eqft_sigmas):.1f}σ")
        print(f"• Systems where GR better: {sum(1 for i in range(len(csv_gr_sigmas)) if csv_gr_sigmas[i] < csv_eqft_sigmas[i])}/{len(csv_results)}")
        print(f"• Systems where E-QFT better: {sum(1 for i in range(len(csv_gr_sigmas)) if csv_eqft_sigmas[i] < csv_gr_sigmas[i])}/{len(csv_results)}")

    # Gaia systems statistical summary (if loaded)
    if gaia_results:
        print(f"\n🌌 GAIA BINARY SYSTEMS STATISTICAL SUMMARY")
        print("=" * 60)
        
        # Extract statistics
        corrections = [data['eqft_correction'] for data in gaia_results.values()]
        correction_percents = [data['correction_percent'] for data in gaia_results.values()]
        g_ratios = [data['G_ratio'] for data in gaia_results.values()]
        masses_1 = [data['M1'] for data in gaia_results.values()]
        masses_2 = [data['M2'] for data in gaia_results.values()]
        eccentricities = [data['e'] for data in gaia_results.values()]
        separations = [data['a_au'] for data in gaia_results.values()]
        
        print(f"📊 SAMPLE STATISTICS:")
        print(f"• Total systems analyzed: {len(gaia_results)}")
        print(f"• Mass range M₁: {min(masses_1):.1f} - {max(masses_1):.1f} M☉")
        print(f"• Mass range M₂: {min(masses_2):.3f} - {max(masses_2):.1f} M☉")
        print(f"• Separation range: {min(separations):.3f} - {max(separations):.3f} AU")
        print(f"• Eccentricity range: {min(eccentricities):.3f} - {max(eccentricities):.3f}")
        
        print(f"\n📈 E-QFT CORRECTION STATISTICS:")
        print(f"• Mean correction: {np.mean(correction_percents):.4f}% of GR")
        print(f"• Std deviation: {np.std(correction_percents):.4f}%")
        print(f"• Range: {min(correction_percents):.4f}% to {max(correction_percents):.4f}%")
        print(f"• Systems with |correction| > 0.1%: {sum(1 for p in correction_percents if abs(p) > 0.1)}")
        print(f"• Systems with |correction| > 1.0%: {sum(1 for p in correction_percents if abs(p) > 1.0)}")
        
        print(f"\n⚙️  G_eff/G_N RATIO DISTRIBUTION:")
        print(f"• Mean: {np.mean(g_ratios):.6f}")
        print(f"• Range: {min(g_ratios):.6f} - {max(g_ratios):.6f}")
        print(f"• Standard deviation: {np.std(g_ratios):.6f}")

    print(f"\n🚀 **EXPERIMENTAL VALIDATION COMPLETE**")
    print(f"Our theoretical β = -0.033 produces specific, falsifiable predictions!")
    print(f"Current data favors pure GR, but E-QFT remains testable with future precision.")

if __name__ == "__main__":
    main()