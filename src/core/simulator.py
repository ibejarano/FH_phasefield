import time
import os
import logging
import numpy as np
from mpi4py import MPI
from .model import HydraulicFractureModel
from .config import SimulationConfig
from src.output_utils import write_output, create_xml_output

logger = logging.getLogger(__name__)

class HydraulicSimulator:
    def __init__(self, model: HydraulicFractureModel, config: SimulationConfig):
        self.model = model
        self.config = config
        
        # Prepare output directory
        if not os.path.exists(config.case_dir):
            os.makedirs(config.case_dir)
            
        self.out_xml = create_xml_output(config.case_dir)
        self.csv_path = os.path.join(config.case_dir, "results.csv")
        
    def run(self, protocol=None):
        """
        Main simulation execution loop.
        Args:
            protocol: Optional function protocol(t) -> Q(t). If None, uses config.Q0.
        """
        mode = "AXISYMMETRIC" if self.config.axisymmetric else ("SYMMETRIC" if self.config.symmetric else "FULL")
        logger.info(f"Starting simulation: {self.config.case_dir} | Mode: {mode}")
        
        t = 0.0
        step = 0
        retry_count = 0
        
        with open(self.csv_path, 'w') as csv_file:
            csv_file.write("time,pressure,volume,wplus,wminus\n")
            
        
        current_dt = self.config.dt
        
        with open(self.csv_path, 'w') as csv_file:
            csv_file.write("time,pressure,volume,wplus,wminus\n")
            
            while t < self.config.t_max:
                step += 1
                
                if self.config.adaptive_time:
                     current_dt = min(current_dt * self.config.dt_growth, self.config.dt_max)
                
                # Determine flow rate Q(t)
                t_target = t + current_dt
                
                Q_t = protocol(t_target) if protocol else self.config.Q0
                
                dV = current_dt * Q_t
                if self.config.symmetric:
                    dV /= 2.0
                
                try:
                    pressure, volume = self.model.solve_time_step(dV)
                    
                    # Adaptive Check
                    increase_dt = True
                    if self.config.adaptive_time:
                        phi_new = self.model.phase.get().vector().get_local()
                        phi_old = self.model.phase.get_old().vector().get_local()
                        dphi = np.max(np.abs(phi_new - phi_old))
                        
                        if dphi > 0.5: # Hard limit
                             raise RuntimeError(f"Step REJECTED: Large phase jump {dphi:.2e}")
                             
                        if dphi > self.config.dphi_max:
                             # If step is too large, reject and shrink
                             # But only if we can shrink further
                             if current_dt > self.config.dt_min:
                                 raise RuntimeError(f"Step REJECTED: dPhi {dphi:.2e} > {self.config.dphi_max}")
                             else:
                                 logger.warning("dPhi limit exceeded but dt is at minimum. Accepting step.")
                        
                    # Commit
                    self.model.commit_state()
                    
                    # Success branch
                    t = t_target
                    retry_count = 0
                    
                    # Compute fracture opening
                    w_plus, w_minus = self.model.compute_openings()
                    
                    # Log progress to CSV
                    csv_file.write(f"{t},{pressure},{volume},{w_plus},{w_minus}\n")
                    
                    # Periodic output
                    if step % self.config.output_freq == 0:
                        self.save_results(t)
                        
                    if step % self.config.store_freq == 0:
                        csv_file.flush()
                    
                    if MPI.COMM_WORLD.rank == 0:
                        dphi_str = f"| dPhi {dphi:.2e}" if self.config.adaptive_time else ""
                        logger.info(f"Step {step:4d} | Time {t:.2e}s | dt {current_dt:.2e} {dphi_str} | P {pressure:8.2f} | V {volume:.2e}")
                
                except Exception as e:
                    if self.config.adaptive_time:
                        logger.warning(f"Step {step} retry due to: {e}")
                        retry_count += 1
                        
                        if retry_count > 5:
                            logger.error("FATAL: Failed trying to find the delta time (too many retries).")
                            self.save_results(t)
                            break
                            
                        current_dt *= self.config.dt_shrink
                        if current_dt < self.config.dt_min:
                            logger.error("FATAL: Minimum time step reached.")
                            self.save_results(t)
                            break
                        
                        step -= 1
                        continue
                    else:
                        logger.error(f"FATAL: Simulation failed at step {step}: {e}")
                        self.save_results(t)
                        break

                
        logger.info("Simulation complete.")

    def save_results(self, t: float):
        """Writes current fields to XML/VTU format."""
        state = self.model.get_current_state()
        write_output(
            self.out_xml, 
            state['u'], 
            state['phi'], 
            state['stress'], 
            t
        )
