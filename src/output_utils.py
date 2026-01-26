from dolfin import File, TimeSeries, XDMFFile, project, conditional, gt, Constant

def create_output_files(caseDir):
    xdmf = XDMFFile(f"{caseDir}/output.xdmf")
    xdmf.parameters["flush_output"] = True
    xdmf.parameters["functions_share_mesh"] = True

    u_ts = TimeSeries(f"{caseDir}/u_series")
    phi_ts = TimeSeries(f"{caseDir}/phi_series")

    return xdmf, u_ts, phi_ts

def write_output(xdmf, u, phi, stress, t):
    # Mask displacement where fracture exists (phi ~ 1) to avoid large visual artifacts
    # We use 0.9 as the cutoff for visualization purposes
    u_masked = project(conditional(gt(phi, 0.9), Constant((0.0, 0.0)), u), u.function_space())
    u_masked.rename("displacement", "displacement")
    
    xdmf.write(u_masked, t)
    xdmf.write(phi, t)
    xdmf.write(stress, t)

def store_time_series(u_ts, phi_ts, u, phi, t):
    u_ts.store(u.vector(), t)
    phi_ts.store(phi.vector(), t)

def create_xml_output(case_dir):
    xdmf = XDMFFile(f"{case_dir}/output.xdmf")
    xdmf.parameters["flush_output"] = True
    xdmf.parameters["functions_share_mesh"] = True
    return xdmf