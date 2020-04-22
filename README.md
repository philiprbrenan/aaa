<h1>Table of contents</h1>
<p><table cellpadding=10>

<tr><td><a href="#01f871a6ce267429760842dc197ac7ca">Eich(HeatLoad)</a><td>

    Class for evaluating and storing the heat load at the up-stream and down-stream

    positions using an Eich-like single exponential fit.



    Member functions:



        calculate_heat_flux_density: Function for calculating the heat-flux given

        a set of input parameters.

<tr><td><a href="#22562053c67373bd208c954c298ee7c3">FieldLine()</a><td>

    Class for following a magnetic field line given a Fiesta equilibrium.



    Member functions:

        follow_field_in_plane(p_0, max_length, max_points, rtol).

<tr><td><a href="#0513469fc7c86890e4049e46a7ab6fb8">Fiesta()</a><td>

    Class for tracing the magnetic field lines given a FIESTA equlibrium.

    This uses a limited set of FIESTA data, making it backward compatible with older FIESTA files,

    but does not have CHERAB functionallity



    :param str filename: the path to the FIESTA MATLAB save file.



    :ivar VectorFunction3D b_field: A 3D vector function of the magnetic field.



    member functions:.

<tr><td><a href="#0e1c0cfbe53c997a5add13b36f7bbe05">HESELdata()</a><td>

    Class for loading and storing the data output from HESEL stored in an hdf5 file.

    All fields are in SI units unless otherwise specified.



    Member functions:

        _load_file()

        _close_file()

        _evaluate_electron_advection_field()

        _evaluate_ion_advection_field()

        _evaluate_electron_conduction_field()

        _evaluate_ion_conduction_field()

        evaluate_parallel_heat_fluxes()

        get_lcfs_values()

        calculate_lambda_q(case='')

        get_probe_positions()

        _load_probe_data(field_name)

        get_profiles_from_probes()

        load_2d_animation_fields()

        animate_2d_field(fieldname='', show_animation=False, save_movie=True)

        animate_1d_field(fieldname='', show_animation=False, save_movie=True).

<tr><td><a href="#92c15be9f555c5f3e8d6c3bf8296bb44">HESELparams()</a><td>

    Class for loading parameters from the HESEL .hdf5 file

    All parameters are in SI units unless otherwise specified.



    Member functions:

        read_hesel_output(self, file)

        calculate_sound_speed(self)

        calculate_ion_gyrofrequency()

        calculate_ion_gyroradius()

        calculate_debye_length()

        calculate_coulomb_logarithm()

        get_x_axis_rhos()

        get_y_axis_rhos()

        get_lcfs_index()

        get_wall_index()

        get_x_axis_probes().

<tr><td><a href="#1b0b1d0a50943a7bf227b2721b9ae77f">HeatLoad()</a><td>

    Base function for storing heat flux profiles.



    Member functions:

        __call__(s) : integrate the heat-flux over s



        calculate_heat_flux_density

        set_coordinates(s_in)

        set_edge_radius(radius_in)

        get_local_coordinates()

        get_global_coordinates()

        calculate_heat_power()

        plot_heat_power_density().

<tr><td><a href="#8ecd2f7bc5142547fcef6ceccbddd4a6">InterfaceSurface</a><td>

    Class for mapping power from an interface surface in the divertor onto the walls.



    :param Point2D point_a: A 2D point representing the start of the interface surface.

    :param Point2D point_b: A 2D point representing the end of the interface surface.

    :param ndarray power_profile: A an array of power values representing the power

      profile along the interface surface. These points are equally spaced and will be

      re-normalised to give an integral of one over the surface.

<tr><td><a href="#4ab17cc669860a18b5efd64996345288">ParticlePath()</a><td>

    Class for determining the path of a charged particle in a magnetic field

    given a Fiesta equilibrium.



    Member functions:

        follow_path(p_0, max_time, max_points, rtol, break_at_limiter).

<tr><td><a href="#2a84ce044b53bb3902f080782d296b09">Vector2()</a><td>Vectors in 2 dimensional euclidean space.

<tr><td><a href="#81c15fb95b1c3422ceb5b2b85c7f9668">add_resource.py</a><td>

<tr><td><a href="#afea94043c155f8b6a9c8265b54802f8">calculate_angle.py</a><td>

<tr><td><a href="#2bc2db0b2b880fd2ffc72ac6bdb18365">equdsk.py</a><td>

<tr><td><a href="#f7f0f2db6aa5ccab3d80ff8536dbf6f7">field_line_projection.py</a><td>

<tr><td><a href="#431e5d48b3ab6a117ee828bfec778106">getOption.py</a><td>

<tr><td><a href="#9386ee3547cc05f3bdb7852a4e7fbdfa">intersections.py</a><td>

<tr><td><a href="#b53b675887a97eff5d9fe276b5c074bb">load_pickle.py</a><td>

<tr><td><a href="#e02daeaf1215d7e93399246760393699">machine_configuration.py</a><td>

<tr><td><a href="#225106ddf5f55966717a2e58dbd07ea3">map_field_lines.py</a><td>

<tr><td><a href="#af5a643bca05127c7079d216ab5fd94c">midplane_power.py</a><td>

<tr><td><a href="#56845d841c468ee501c30899177b17c6">midplane_to_divertor_power.py</a><td>

<tr><td><a href="#9fb2a57520e98a6c2e8897b9609268d9">parser.py</a><td>

<tr><td><a href="#c49c1165dd9fa01598d89310eece1054">physics_constants.py</a><td>

<tr><td><a href="#cc87324b013b0c555816970d52f4e4d4">project_heat_flux.py</a><td>

<tr><td><a href="#537ca69b729227582c80ebad18406ae1">psi_map_projection.py</a><td>

<tr><td><a href="#4763e95919953c2a5b089bdd3fafb4d5">remove_resource.py</a><td>

<tr><td><a href="#7db598b19ebf90a0f9549002f7f2d962">resource_manager.py</a><td>

<tr><td><a href="#6612d973c9ef53d7152af2c08956fc7a">update_resource.py</a><td>

<tr><td><a href="#9ebf128f7c57b374928a8d28e3bc06f4">vitarunner.py</a><td>

</table>

<h1 id="01f871a6ce267429760842dc197ac7ca">Class: Eich(HeatLoad)</h1>
<p>

    Class for evaluating and storing the heat load at the up-stream and down-stream

    positions using an Eich-like single exponential fit.



    Member functions:



        calculate_heat_flux_density: Function for calculating the heat-flux given

        a set of input parameters.</p>

<h2>__init__</h2>


        Inputs for the class are



        <p><b>Parameters</b>:

        lambda_q : float, optional

            The scrape-off layer power fall-off length. The default is 0.0015 m.

        S : float, optional

            Spreading factor. The default is 0.

        r0_lfs : float, optional

            The position of the low field side LCFS. The default is 0.

        r0_hfs : float, optional

            The position of the high field side LCFS. The default is 0.



        <p><b>Returns</b>:

        None.

.
<h4>Parameters</h4>
<table cellpadding=10>
<tr><th>Name<th>Description
<tr><td><b><i>self</i></b><td>lambda_q=1.5e-3

</table>

<h4>Examples</h4>
<pre language="python">

</pre>
<hr>
<h2>calculate_heat_flux_density</h2>


        Function for evaluating the heat flux density using an Eich-like

        exponential at the mid-plane and an Eich function at the divertor.



        The heat flux for the different positions is given by:



            lfs-mp:

                lambda_q = q_0*exp[-s/lambda_q] for s > 0



            hfs-mp:

                lambda_q = q_0*exp[-s/lambda_q]/f_x_in_out*r0_lfs/r0_hfs

                                                for s > 0 and s < d_rsep

                lambda_q = 0                    for s > d_rsep



            lfs:

                lambda_q = q_0/2*exp[S^2/(2*lambda_q*f_x)^2 - s/(lambda_q*f_x)]

                                *erfc(S/(2*lambda_q*f_x) - s/S)

            hfs:

                So far, this has not been implemented properly, so it simply

                gives the same output as hfs-mp



        <p><b>Parameters</b>:

        where : string, optional

            A string with the position to evaluate.

            The options are:

                "lfs-mp"

                "hfs-mp"

                "lfs"

                "hfs"

            The default is "lfs-mp".



        <p><b>Returns</b>:

        None.

.
<h4>Parameters</h4>
<table cellpadding=10>
<tr><th>Name<th>Description
<tr><td><b><i>self</i></b><td>where="lfs-mp"

</table>

<h4>Examples</h4>
<pre language="python">

</pre>
<hr>
<h1 id="22562053c67373bd208c954c298ee7c3">Class: FieldLine()</h1>
<p>

    Class for following a magnetic field line given a Fiesta equilibrium.



    Member functions:

        follow_field_in_plane(p_0, max_length, max_points, rtol).</p>

<h2>__init__</h2>

<h4>Parameters</h4>
<table cellpadding=10>
<tr><th>Name<th>Description
<tr><td><b><i>self</i></b><td>fiesta

</table>

<h4>Examples</h4>
<pre language="python">

</pre>
<hr>
<h2>dx_dl</h2>


            The function describing the ode to solve in order to track the magnetic field lines



            <p><b>Input</b>: l_dist, np.array with the distance along the magnetic field-line

                   x_vec,  vector with the R, phi and Z initial positions



            <p><b>Return</b>: dx_dl_rhs, the right-hand side of the ode to solve.
<h4>Parameters</h4>
<table cellpadding=10>
<tr><th>Name<th>Description
<tr><td><b><i>_l_dist</i></b><td>x_vec

</table>

<h4>Examples</h4>
<pre language="python">

</pre>
<hr>
<h2>event</h2>


            Function for determining whether the solution to the ODE passes a wall,

            which terminates the ODE solver



            <p><b>Input</b>: l_dist, np.array with the distance along the magnetic field-line

                   x_vec,  vector with the R, phi and Z initial positions



            <p><b>Return</b>: intersect_wall, returns 0 if any wall surface is intersected,

                                    otherwise returns a float (the event function of

                                    solve_ivp only looks for events = 0).
<h4>Parameters</h4>
<table cellpadding=10>
<tr><th>Name<th>Description
<tr><td><b><i>_l_dist</i></b><td>x_vec

</table>

<h4>Examples</h4>
<pre language="python">

</pre>
<hr>
<h1 id="0513469fc7c86890e4049e46a7ab6fb8">Class: Fiesta()</h1>
<p>

    Class for tracing the magnetic field lines given a FIESTA equlibrium.

    This uses a limited set of FIESTA data, making it backward compatible with older FIESTA files,

    but does not have CHERAB functionallity



    :param str filename: the path to the FIESTA MATLAB save file.



    :ivar VectorFunction3D b_field: A 3D vector function of the magnetic field.



    member functions:.</p>

<h2>__init__</h2>

<h4>Parameters</h4>
<table cellpadding=10>
<tr><th>Name<th>Description
<tr><td><b><i>self</i></b><td>filename

<tr><td><b><i>self</i></b><td>filename

</table>

<h4>Examples</h4>
<pre language="python">

</pre>
<hr>
<h2>b_field</h2>


        Function for getting the magnetic a 3D vector function of the magnetic 

        field.
<h4>Parameters</h4>
<table cellpadding=10>
<tr><th>Name<th>Description
<tr><td><b><i>self</i></b><td>

<tr><td><b><i>self</i></b><td>

</table>

<h4>Examples</h4>
<pre language="python">

</pre>
<hr>
<h2>get_midplane_lcfs</h2>


        Function for getting the inner and outer radial position of the LCFS at the midplane



        <p><b>Input</b>: self,  a reference to the object itself

               psi_p, the flux surface of the LCFS, standard is psi_p = 1.005 (otherwise the field-line

                      is located inside the LCFS)



        <p><b>Return</b>: Rcross, a list with the outer and inner radial position of the mid-plane LCFS.
<h4>Parameters</h4>
<table cellpadding=10>
<tr><th>Name<th>Description
<tr><td><b><i>self</i></b><td>psi_p=1.005

<tr><td><b><i>self</i></b><td>psi_p=1.005

</table>

<h4>Examples</h4>
<pre language="python">

</pre>
<hr>
<h2>read_fiesta_model</h2>


        Function for reading the FIESTA equilibrium data from a .mat file



        <p><b>Input</b>: self, a reference the object itself



        <p><b>Output</b>: self.r_limiter, a numpy array with the radial coordinates of the vessel limits

                self.z_limiter, a numpy array with the vertical coordinates of the vessel limits

                self.r_vec,     a numpy array with the radial grid coordinates

                self.z_vec,     a numpy array with the vertical grid coordinates

                self.psi_n,

                self.b_r,       a numpy array with the radial magnetic field component

                self.b_z        a numpy array with the vertical field component

                self.b_phi,     a numpy array with the toroidal field component

                self.b_theta,   a numpy array with the poloidal magnetic field component

                self.i_rod      a float with the current in the rod.
<h4>Parameters</h4>
<table cellpadding=10>
<tr><th>Name<th>Description
<tr><td><b><i>self</i></b><td>

<tr><td><b><i>self</i></b><td>

</table>

<h4>Examples</h4>
<pre language="python">

</pre>
<hr>
<h2>to_cherab_equilibrium</h2>


        Function for converting this Fiesta object to a CHERAB equilibrium.



        rtype: EFITEquilibrium.
<h4>Parameters</h4>
<table cellpadding=10>
<tr><th>Name<th>Description
<tr><td><b><i>self</i></b><td>

</table>

<h4>Examples</h4>
<pre language="python">

</pre>
<hr>
<h1 id="0e1c0cfbe53c997a5add13b36f7bbe05">Class: HESELdata()</h1>
<p>

    Class for loading and storing the data output from HESEL stored in an hdf5 file.

    All fields are in SI units unless otherwise specified.



    Member functions:

        _load_file()

        _close_file()

        _evaluate_electron_advection_field()

        _evaluate_ion_advection_field()

        _evaluate_electron_conduction_field()

        _evaluate_ion_conduction_field()

        evaluate_parallel_heat_fluxes()

        get_lcfs_values()

        calculate_lambda_q(case='')

        get_probe_positions()

        _load_probe_data(field_name)

        get_profiles_from_probes()

        load_2d_animation_fields()

        animate_2d_field(fieldname='', show_animation=False, save_movie=True)

        animate_1d_field(fieldname='', show_animation=False, save_movie=True).</p>

<h2>__del__</h2>

<h4>Parameters</h4>
<table cellpadding=10>
<tr><th>Name<th>Description
<tr><td><b><i>self</i></b><td>

</table>

<h4>Examples</h4>
<pre language="python">

</pre>
<hr>
<h2>__init__</h2>


        Initialise the HESEL data object as None-types



        <p><b>Input</b>: filename, the name of the file to be loaded

               ratio,    the ratio of timesteps to be filtered out to make sure

                         a turbulent steady-state has been reached,

                         e.g., ratio = 0.25 means the first 25% of the data

                         is discarded.
<h4>Parameters</h4>
<table cellpadding=10>
<tr><th>Name<th>Description
<tr><td><b><i>self</i></b><td>filename

</table>

<h4>Examples</h4>
<pre language="python">

</pre>
<hr>
<h2>_close_file</h2>


        If the file is open, close the file.
<h4>Parameters</h4>
<table cellpadding=10>
<tr><th>Name<th>Description
<tr><td><b><i>self</i></b><td>

</table>

<h4>Examples</h4>
<pre language="python">

</pre>
<hr>
<h2>_evaluate_electron_advection_field</h2>


        Function for evaluating the electron advection for each point in time and

        space:



            pe_adv = 3/2*1/tau_d*sqrt((T_i + T_e)/(T_i0 + T_e0))*p_e,



        where the sqrt((T_i + T_e)/(T_i0 + T_e0))*p_e is to take local variations

        into account and 3/2 is from the normalisation used in HESEL.



        <p><b>Input</b>: self,



        <p><b>Return</b>: electron_advecton_mw, a numpy array with the ion advection term for

                each point in space and time [MW].
<h4>Parameters</h4>
<table cellpadding=10>
<tr><th>Name<th>Description
<tr><td><b><i>self</i></b><td>

</table>

<h4>Examples</h4>
<pre language="python">

</pre>
<hr>
<h2>_evaluate_electron_conduction_field</h2>


        Function for evaluating the electron conduction for each point in time and

        space:



            pe_cond = 3/2*1/tau_{SH,e}*(T_e/T_e0)^(5/2)*p_e,



        where the (T_e/T_e0)^(5/2)*p_e is to take local variations into account

        and 3/2 is from the normalisation used in HESEL.



        <p><b>Input</b>: self,



        <p><b>Return</b>: electron_conduction_mw, a numpy array with the electron conduction

                term for each point in space and time [MW].
<h4>Parameters</h4>
<table cellpadding=10>
<tr><th>Name<th>Description
<tr><td><b><i>self</i></b><td>

</table>

<h4>Examples</h4>
<pre language="python">

</pre>
<hr>
<h2>_evaluate_ion_advection_field</h2>


        Function for evaluating the ion advection for each point in time and

        space:



            pi_adv = 3/2*1/tau_d*sqrt((T_i + T_e)/(T_i0 + T_e0))*p_i,



        where the sqrt((T_i + T_e)/(T_i0 + T_e0))*p_i is to take local variations

        into account and 3/2 is from the normalisation used in HESEL.



        <p><b>Input</b>: self,



        <p><b>Return</b>: ion_advecton_mw, a numpy array with the ion advection term for

                each point in space and time [MW].
<h4>Parameters</h4>
<table cellpadding=10>
<tr><th>Name<th>Description
<tr><td><b><i>self</i></b><td>

</table>

<h4>Examples</h4>
<pre language="python">

</pre>
<hr>
<h2>_evaluate_ion_conduction_field</h2>


        Function for evaluating the electron conduction for each point in time and

        space:



            pi_cond = 3/2*1/tau_{SH,i}*(T_i/T_e0)^(5/2)*p_i,



        where the (T_i/T_e0)^(5/2)*p_i is to take local variations into account

        (T_i is normalised to T_e0 in HESEL) and 3/2 is from the normalisation used in HESEL.



        <p><b>Input</b>: self,



        <p><b>Return</b>: electron_conduction_mw, a numpy array with the electron conduction

                term for each point in space and time [MW].
<h4>Parameters</h4>
<table cellpadding=10>
<tr><th>Name<th>Description
<tr><td><b><i>self</i></b><td>

</table>

<h4>Examples</h4>
<pre language="python">

</pre>
<hr>
<h2>_evaluate_parallel_heat_flux_electron_advection</h2>


        Function for evaluate the contribution to the parallel heat flux from

        the ion advection:



            q_parallel_e_adv = a*<p_e/tau_d>_t,



        where a is the device minor radius, tau_d is the parallel advection loss term

        and <>_t denotes a temporal average. There are fewer probe radial points than

        1D radial points in HESEL, so the heat flux is interpolated after it is evaluated.



        <p><b>Input</b>: self,



        <p><b>Return</b>: q_parallel_e_adv_mw, a numpy array with the electron advection part

                                     of the parallel heat flux.
<h4>Parameters</h4>
<table cellpadding=10>
<tr><th>Name<th>Description
<tr><td><b><i>self</i></b><td>

</table>

<h4>Examples</h4>
<pre language="python">

</pre>
<hr>
<h2>_evaluate_parallel_heat_flux_electron_conduction</h2>


        Function for evaluate the contribution to the parallel heat flux from

        the ion advection:



            q_parallel_e_cond = a*<p_e/tau_{SH,e}>_t,



        where a is the device minor radius, tau_{SH,e} is the parallel conduction loss term

        and <>_t denotes a temporal average. There are fewer probe radial points than

        1D radial points in HESEL, so the heat flux is interpolated after it is evaluated.



        <p><b>Input</b>: self,



        <p><b>Return</b>: q_parallel_e_cond_mw, a numpy array with the electron conduction part

                                      of the parallel heat flux.
<h4>Parameters</h4>
<table cellpadding=10>
<tr><th>Name<th>Description
<tr><td><b><i>self</i></b><td>

</table>

<h4>Examples</h4>
<pre language="python">

</pre>
<hr>
<h2>_evaluate_parallel_heat_flux_ion_advection</h2>


        Function for evaluate the contribution to the parallel heat flux from

        the ion advection:



            q_parallel_i_adv = a*<p_i/tau_d>_t,



        where a is the device minor radius, tau_d is the parallel advection loss term

        and <>_t denotes a temporal average. There are fewer probe radial points than

        1D radial points in HESEL, so the heat flux is interpolated after it is evaluated.



        <p><b>Input</b>: self,



        <p><b>Return</b>: q_parallel_i_adv_mw, a numpy array with the ion advection part

                                     of the parallel heat flux.
<h4>Parameters</h4>
<table cellpadding=10>
<tr><th>Name<th>Description
<tr><td><b><i>self</i></b><td>

</table>

<h4>Examples</h4>
<pre language="python">

</pre>
<hr>
<h2>_evaluate_parallel_heat_flux_ion_conduction</h2>


        Function for evaluate the contribution to the parallel heat flux from

        the ion advection:



            q_parallel_i_cond = a*<p_i/tau_{SH,i}>_t,



        where a is the device minor radius, tau_{SH,i} is the parallel conduction loss term

        and <>_t denotes a temporal average. There are fewer probe radial points than

        1D radial points in HESEL, so the heat flux is interpolated after it is evaluated.



        <p><b>Input</b>: self,



        <p><b>Return</b>: q_parallel_i_cond_mw, a numpy array with the ion conduction part

                                      of the parallel heat flux.
<h4>Parameters</h4>
<table cellpadding=10>
<tr><th>Name<th>Description
<tr><td><b><i>self</i></b><td>

</table>

<h4>Examples</h4>
<pre language="python">

</pre>
<hr>
<h2>_load_file</h2>


        Load the file and the HESEL parameters from HESELparameters.py.
<h4>Parameters</h4>
<table cellpadding=10>
<tr><th>Name<th>Description
<tr><td><b><i>self</i></b><td>

</table>

<h4>Examples</h4>
<pre language="python">

</pre>
<hr>
<h2>_load_probe_data</h2>


        Function for loading the probe synthetic HESEL probe data at the radial position

        radial_probe_position = 0.0 and at poloidal positions that are in increments of

        5 rho_s from poloidal_probe_position = 0.0, which is where we assume that the

        synthetic data is mutually independent.



        The data is sliced so the first self.ratio*n_t time-points are filtered away from the

        signals of each synthetic probe. The output numpy array then consists of the

        concatenated synthetic data from the probes at different poloidal positions as

        specified above.



        <p><b>Input</b>: self,

               field_name,  a string with the name of the field to load. The options are

                            'temperature' for the electron temperature,

                            'temperature_i' for the ion temperature, and

                            'density' for the plasma density (the plasma is assumed quasi-

                            neutral)



        <p><b>Return</b>: field_data, a numpy array with the specified field in SI units.
<h4>Parameters</h4>
<table cellpadding=10>
<tr><th>Name<th>Description
<tr><td><b><i>self</i></b><td>field_name

</table>

<h4>Examples</h4>
<pre language="python">

</pre>
<hr>
<h2>animate_1d_field</h2>


        Function for creating videos of the 1D fields



        <p><b>Input</b>: self

               fieldname, a string with the name of the field to make the video of

                          the options are: 'n', 'Te', 'Ti', 'Pe', 'Pi'

                          default is '', which returns an error

               show_animation, a boolean determining whether or not to show the

                               animation as it is being created

                               default is False

               save_movie, a boolean for whether or not the movie should be saved,

                           default is True



        <p><b>Output</b>: An animation with the specified field,

                saved in the working directory if save_movie == True.
<h4>Parameters</h4>
<table cellpadding=10>
<tr><th>Name<th>Description
<tr><td><b><i>self</i></b><td>fieldname=''

</table>

<h4>Examples</h4>
<pre language="python">

</pre>
<hr>
<h2>animate_2d_field</h2>


        Function for creating videos of the 2D fields



        <p><b>Input</b>: self,

               fieldname, a string with the name of the field to make the video of

                          the options are: 'n', 'Te', 'Ti', 'Pe', 'Pi', 'phi', 'omega'

                          default is '', which returns an error

               show_animation, a boolean determining whether or not to show the

                               animation as it is being created

                               default is False

               save_movie, a boolean for whether or not the movie should be saved,

                           default is True



        <p><b>Output</b>: An animation with the specified field,

                saved in the working directory if save_movie == True.
<h4>Parameters</h4>
<table cellpadding=10>
<tr><th>Name<th>Description
<tr><td><b><i>self</i></b><td>fieldname=''

</table>

<h4>Examples</h4>
<pre language="python">

</pre>
<hr>
<h2>calculate_lambda_q</h2>


        Function for calculating lambda_q as the weighted average position:



            lambda_q = (int_0^l_x x*q_parallel(x) dx)/(int_0^l_x q_parallel(x) dx),



        where l_x is the length of the domain, x is the radial coordinates and

        q_parallel(x) is the parallel heat flux profile



        <p><b>Input</b>: self,

               case, a string with which lambda_q to evaluate. Can be:

                     'q_tot', only evaluate lambda_q on q_parallel_tot

                     'q_adv_e', only evaluate lambda_q on q_parallel_e_adv

                     'q_con_e', only evaluate lambda_q on q_parallel_e_con

                     'q_adv_i', only evaluate lambda_q on q_parallel_i_adv

                     'q_con_i', only evaluate lambda_q on q_parallel_i_con

                     If nothing is stated, all of them are evaluated



        <p><b>Output</b>: self.lambda_q_tot,   a float with lambda_q_tot if 'q_tot' is specified

                self.lambda_q_e_adv, a float with lambda_q_e_adv if 'q_adv_e' is specified

                self.lambda_q_e_con, a float with lambda_q_e_con if 'q_con_e' is specified

                self.lambda_q_i_adv, a float with lambda_q_i_adv if 'q_adv_i' is specified

                self.lambda_q_i_con, a float with lambda_q_i_con if 'q_con_i' is specified



                All of the above are output if anything else or nothing is specified.
<h4>Parameters</h4>
<table cellpadding=10>
<tr><th>Name<th>Description
<tr><td><b><i>self</i></b><td>case=''

</table>

<h4>Examples</h4>
<pre language="python">

</pre>
<hr>
<h2>evaluate_parallel_heat_fluxes</h2>


        Function for evaluating the parallel heat fluxes



        <p><b>Input</b>: self

        <p><b>Output</b>: self.q_parallel_tot,   the total parallel heat flux profile in MW/m^2

                self.q_parallel_e_con, the parallel heat flux profile from the

                                       electron conduction in MW/m^2

                self.q_parallel_e_adv, the parallel heat flux profile from the

                                       electron advection in MW/m^2

                self.q_parallel_i_con, the parallel heat flux profile from the

                                       ion conduction in MW/m^2

                self.q_parallel_i_adv, the parallel heat flux profile from the

                                       ion advection in MW/m^2.
<h4>Parameters</h4>
<table cellpadding=10>
<tr><th>Name<th>Description
<tr><td><b><i>self</i></b><td>

</table>

<h4>Examples</h4>
<pre language="python">

</pre>
<hr>
<h2>get_lcfs_values</h2>


        Function for getting calculating the average LCFS values of the plasma parameters.

        The 1D fields in HESEL have better spatial resolution than the synthetic probes,

        so we use those to evaluate the values at the LCFS



        <p><b>Input</b>:  self,



        <p><b>Output</b>: self.n_lcfs,       a float with the average density at the LCFS

                self.te_lcfs_ev,   a float with the the electron temperature at the LCFS in eV

                self.ti_lcfs_ev,   a float with the the ion temperature at the LCFS in eV

                self.grad_pe_lcfs, a float with the the gradient of the electron pressure

                                   at the LCFS

                self.grad_pi_lcfs, a float with the  the gradient of the ion pressure at the LCFS.
<h4>Parameters</h4>
<table cellpadding=10>
<tr><th>Name<th>Description
<tr><td><b><i>self</i></b><td>

</table>

<h4>Examples</h4>
<pre language="python">

</pre>
<hr>
<h2>get_probe_positions</h2>


        Function for getting the positions of the synthetic probes used in the simulation.

        The information is stored in the file called myprobe.dat, which is converted to

        a string when loading it into python.



        Probe names are specified with @TIP followed by a space '\t', then the radial

        position, a space '\t' and then the poloidal position.

        The end of the segment containing positions ends with a space '\t' and 'hdf5'.



        <p><b>Input</b>:  self,



        <p><b>Output</b>: self.probe_position,  a dictionary with the probename as the key and a list

                                      with the radial and poloidal position of the synthetic probes

                                      [radial_position, poloidal_position].
<h4>Parameters</h4>
<table cellpadding=10>
<tr><th>Name<th>Description
<tr><td><b><i>self</i></b><td>

</table>

<h4>Examples</h4>
<pre language="python">

</pre>
<hr>
<h2>get_profiles_from_probes</h2>


        Function for getting n, t_e and t_i profiles from the synthetic probe

        diagnostic in the HESEL hdf5-file. The synthetic probes have lower spatial

        resolution than the 1D fields, but have a higher temporal resolution.



        <p><b>Input</b>:  self,



        <p><b>Output</b>: self.n_probes, a numpy array with the probe profile data

                               for the density

                self.te_probes, a numpy array with the probe profile data

                                for the electron temperature

                self.ti_probes, a numpy array with the probe profile data

                                for the ion temperature.
<h4>Parameters</h4>
<table cellpadding=10>
<tr><th>Name<th>Description
<tr><td><b><i>self</i></b><td>

</table>

<h4>Examples</h4>
<pre language="python">

</pre>
<hr>
<h2>load_2d_animation_fields</h2>


        Load the 2D fields from the HESEL code from xanimation (the spatial resolution

        is 1/4 that of the full 2D profiles)



        <p><b>Input</b>:  self

        <p><b>Output</b>: self.n_2d, the 2D field of the density in SI units

                self.pe_2d, the 2D field of the electron pressure in SI units

                self.te_2d, the 2D field of the electron temperature in SI units,

                            calculated as P_e/n

                self.pi_2d, the 2D field of the ion pressure in SI units

                self.ti_2d, the 2D field of the ion temperature in SI units,

                            calculated as P_i/n

                self.phi_2d, the 2D field of the potential in SI units.
<h4>Parameters</h4>
<table cellpadding=10>
<tr><th>Name<th>Description
<tr><td><b><i>self</i></b><td>

</table>

<h4>Examples</h4>
<pre language="python">

</pre>
<hr>
<h1 id="92c15be9f555c5f3e8d6c3bf8296bb44">Class: HESELparams()</h1>
<p>

    Class for loading parameters from the HESEL .hdf5 file

    All parameters are in SI units unless otherwise specified.



    Member functions:

        read_hesel_output(self, file)

        calculate_sound_speed(self)

        calculate_ion_gyrofrequency()

        calculate_ion_gyroradius()

        calculate_debye_length()

        calculate_coulomb_logarithm()

        get_x_axis_rhos()

        get_y_axis_rhos()

        get_lcfs_index()

        get_wall_index()

        get_x_axis_probes().</p>

<h2>__init__</h2>

<h4>Parameters</h4>
<table cellpadding=10>
<tr><th>Name<th>Description
<tr><td><b><i>self</i></b><td>file

</table>

<h4>Examples</h4>
<pre language="python">

</pre>
<hr>
<h2>calculate_coulomb_logarithm</h2>


        Function for calculating the Coulomb logarithm:



            log(Lambda_coulomb) = log(12 pi n_0 lambda_debye^3/Z),



        n0 is the reference plasma density, lambda_debye is the plasma debye length and

        Z is the plasma charge number



        <p><b>Input</b>: self,



        <p><b>Return</b>: log(Lambda_coulomb), the plasma coulomb logarithm.
<h4>Parameters</h4>
<table cellpadding=10>
<tr><th>Name<th>Description
<tr><td><b><i>self</i></b><td>

</table>

<h4>Examples</h4>
<pre language="python">

</pre>
<hr>
<h2>calculate_debye_length</h2>


        Function for calculating the plasma debye length:



            lambda_debye = sqrt(epsilon Te0/(n_0 e)),



        where epsilon is the vacuum permittivity, Te is the reference electron temperature in eV,

        n0 is the reference plasma density and e is the elementary charge



        <p><b>Input</b>: self,



        <p><b>Return</b>: lambda_debye, the plasma debye length.
<h4>Parameters</h4>
<table cellpadding=10>
<tr><th>Name<th>Description
<tr><td><b><i>self</i></b><td>

</table>

<h4>Examples</h4>
<pre language="python">

</pre>
<hr>
<h2>calculate_ion_gyrofrequency</h2>


        Function for calculating the ion gyrofrequency at electron temperature:



            omega_ci = Z e B_0/m_i,



        where Z is the charge number, e is the elementary charge,

        B_0 is the magnetic field strength at the outer midplane and

        m_i is the ion mass.



        <p><b>Input</b>: self,



        <p><b>Return</b>: omega_ci, the ion sound speed at electron temperature.
<h4>Parameters</h4>
<table cellpadding=10>
<tr><th>Name<th>Description
<tr><td><b><i>self</i></b><td>

</table>

<h4>Examples</h4>
<pre language="python">

</pre>
<hr>
<h2>calculate_ion_gyroradius</h2>


        Function for calculating the ion gyrofrequency at electron temperature:



            rho_s = c_s/omega_ci,



        where c_s is the ion sound speed at electron temperature and

        omega_ci is the ion gyrofrequency at electron temperature



        <p><b>Input</b>: self,



        <p><b>Return</b>: rho_s, the ion gyroradius at electron temperature.
<h4>Parameters</h4>
<table cellpadding=10>
<tr><th>Name<th>Description
<tr><td><b><i>self</i></b><td>

</table>

<h4>Examples</h4>
<pre language="python">

</pre>
<hr>
<h2>calculate_sound_speed</h2>


        Function for calculating the ion sound speed at electron temperature:



            c_s = sqrt(e Te/m_i),



        where e is the elementary charge, Te is the reference electron temperature in eV and

        m_i is the ion mass.



        <p><b>Input</b>: self,



        <p><b>Return</b>: c_s, the ion sound speed at electron temperature.
<h4>Parameters</h4>
<table cellpadding=10>
<tr><th>Name<th>Description
<tr><td><b><i>self</i></b><td>

</table>

<h4>Examples</h4>
<pre language="python">

</pre>
<hr>
<h2>get_lcfs_index</h2>


        Function for calculating the index of the LCFS position in the HESEL output



        <p><b>Input</b>: self,



        <p><b>Return</b>: lcfs_index, an integer with the index of the LCFS.
<h4>Parameters</h4>
<table cellpadding=10>
<tr><th>Name<th>Description
<tr><td><b><i>self</i></b><td>

</table>

<h4>Examples</h4>
<pre language="python">

</pre>
<hr>
<h2>get_wall_index</h2>


        Function for calculating the index where the wall region starts



        <p><b>Input</b>: self,



        <p><b>Return</b>: wall_index, an integer with the index of where the wall region starts.
<h4>Parameters</h4>
<table cellpadding=10>
<tr><th>Name<th>Description
<tr><td><b><i>self</i></b><td>

</table>

<h4>Examples</h4>
<pre language="python">

</pre>
<hr>
<h2>get_x_axis_probes</h2>


        Function for getting the x-axis for the synthetic probes in HESEL, normalised

        so the LCFS is at 0



        <p><b>Input</b>: self,



        <p><b>Return</b>: x_axis_probes_rhos, a numpy array with the synthetic probe positions.
<h4>Parameters</h4>
<table cellpadding=10>
<tr><th>Name<th>Description
<tr><td><b><i>self</i></b><td>

</table>

<h4>Examples</h4>
<pre language="python">

</pre>
<hr>
<h2>get_x_axis_rhos</h2>


        Function for getting the x-axis of the HESEL output, normalised so 0 is at the LCFS



        <p><b>Input</b>: self,



        <p><b>Return</b>: x_axis, a numpy array with the radial positions of the grid points [rho_s].
<h4>Parameters</h4>
<table cellpadding=10>
<tr><th>Name<th>Description
<tr><td><b><i>self</i></b><td>

</table>

<h4>Examples</h4>
<pre language="python">

</pre>
<hr>
<h2>get_y_axis_rhos</h2>


        Function for getting the y-axis of the HESEL output



        <p><b>Input</b>: self,



        <p><b>Return</b>: y_axis, a numpy array with the poloidal positions of the grid points [rho_s].
<h4>Parameters</h4>
<table cellpadding=10>
<tr><th>Name<th>Description
<tr><td><b><i>self</i></b><td>

</table>

<h4>Examples</h4>
<pre language="python">

</pre>
<hr>
<h2>read_hesel_output</h2>


        Function for reading the output from the HESEL .hdf5 file.

        All data is in SI units unless otherwise specified.



        <p><b>Input</b>: self,

               file, the loaded HESEL .hdf5 file



        <p><b>Output</b>: self.n_x,              an integer with the number of radial grid points

                self.n_y,              an integer with the number of poloidal grid points

                self.n_t,              an integer with the number of temporal points

                self.outmult,          an integer with the number of timesteps per output

                self.xmin_rhos,        a float with the minimum radial position

                                       of the domain in [rho_s]

                self.xmax_rhos,        a float with the maximum radial position

                                       of the domain in [rho_s]

                self.n_0,              a float with reference plasma density [m^{-3}]

                self.te0_eV,           a float with reference plasma electron temperature [eV]

                self.ti0_eV,           a float with reference plasma ion temperature [eV]

                self.background_n,     a float with the background density [m^{-3}]

                self.background_t,     a float with the background temperatures

                                       (both ele and ion) [eV]

                self.plasma_z,         an integer with the plasma charge number

                self.b0_omp,           a float with the magnitude of the magnetic

                                       field at the OMP [T]

                self.ion_mass_number,  an integer with the ion mass number

                self.minor_radius,     a float with the device minor radius [m]

                self.major_radius,     a float with the device major radius [m]

                self.plasma_q,         a float with the plasma safety factor

                self.parallel_conn_length,      a float with the parallel connection length [m]

                self.parallel_conn_length_wall, a float with the parallel connection length

                                                of the wall region [m]

                self.mach_number,      a float with the parallel mach number at the OMP

                self.edge_width_rhos,  a float with the width of the edge region [rho_s]

                self.sol_width_rhos,   a float with the width of the SOL region [rho_s]

                self.wall_region_width_rhos,    a float with the width of the

                                                wall region [rho_s]

                self.time1d_omega_ci,  a float with the time-step used in HESEL [1/omega_ci]

                self.time2d_omega_ci,  a float with the output time-step used for 2D fields in HESEL

                self.probes_nt,        an integer with the number of output times

                                       for the 1D probe arrays

                self.probes_nx,        an integer with the radial number of probes

                                       in the 1D probe arrays

                self.dx_rhos,          a float with the radial resolution of the domain [rho_s]

                self.dy_rhos,          a float with the poloidal resolution of the domain [rho_s]

                self.adv_p_e_omega_ci, a float with 9/2*1/(omega_ci*tau_s), i.e. the normalised

                                       constant part of the parallel advection term

                self.adv_p_i_omega_ci, a float with 9/2*1/(omega_ci*tau_s), i.e. the normalised

                                       constant part of the parallel advection term

                self.con_p_e_omega_ci, a float with 1/(omega_ci*tau_{SH,e}), i.e. the normalised

                                       const. part of the parallel electron Spitzer-Härm conduction

                self.con_p_i_omega_ci, a float with 1/(omega_ci*tau_{SH,i}), i.e. the const.

                                       normalised part of the parallel ion Spitzer-Härm conduction.
<h4>Parameters</h4>
<table cellpadding=10>
<tr><th>Name<th>Description
<tr><td><b><i>self</i></b><td>file

</table>

<h4>Examples</h4>
<pre language="python">

</pre>
<hr>
<h1 id="1b0b1d0a50943a7bf227b2721b9ae77f">Class: HeatLoad()</h1>
<p>

    Base function for storing heat flux profiles.



    Member functions:

        __call__(s) : integrate the heat-flux over s



        calculate_heat_flux_density

        set_coordinates(s_in)

        set_edge_radius(radius_in)

        get_local_coordinates()

        get_global_coordinates()

        calculate_heat_power()

        plot_heat_power_density().</p>

<h2>__call__</h2>

<h4>Parameters</h4>
<table cellpadding=10>
<tr><th>Name<th>Description
<tr><td><b><i>self</i></b><td>s

</table>

<h4>Examples</h4>
<pre language="python">

</pre>
<hr>
<h2>__init__</h2>

<h4>Parameters</h4>
<table cellpadding=10>
<tr><th>Name<th>Description
<tr><td><b><i>self</i></b><td>r0_lfs=0.

</table>

<h4>Examples</h4>
<pre language="python">

</pre>
<hr>
<h2>calculate_heat_flux_density</h2>


        Function for calculating the heat flux profile. Needs to be implemented

        by the inhering class.



        Raises

        ------

        NotImplementedError



        <p><b>Returns</b>:

        None.

.
<h4>Parameters</h4>
<table cellpadding=10>
<tr><th>Name<th>Description
<tr><td><b><i>self</i></b><td>

</table>

<h4>Examples</h4>
<pre language="python">

</pre>
<hr>
<h2>calculate_heat_power</h2>


        Function for calculating the integral of the power profile as a function

        of s



            heat_power = p_tot/(2*pi*(R_0+a)) = int_0^L_perp q(s) ds



        <p><b>Returns</b>:

        heat_power : float

            The integral of the heat flux profile with respect to s

.
<h4>Parameters</h4>
<table cellpadding=10>
<tr><th>Name<th>Description
<tr><td><b><i>self</i></b><td>

</table>

<h4>Examples</h4>
<pre language="python">

</pre>
<hr>
<h2>get_global_coordinates</h2>


        Function for getting the global coordinates, i.e. the local coordinates

        plus the position of the LCFS for the LFS and the position of the LCFS

        minus the local position on the HFS



        <p><b>Returns</b>:

        _s_global : list or 1-by-n numpy array

            An array with the global coordinates.

.
<h4>Parameters</h4>
<table cellpadding=10>
<tr><th>Name<th>Description
<tr><td><b><i>self</i></b><td>location='lfs'

</table>

<h4>Examples</h4>
<pre language="python">

</pre>
<hr>
<h2>get_local_coordinates</h2>


        Function for getting the local coordinates



        <p><b>Returns</b>:

        _s : list or 1-by-n numpy array

            An array with the local coordinates.

.
<h4>Parameters</h4>
<table cellpadding=10>
<tr><th>Name<th>Description
<tr><td><b><i>self</i></b><td>

</table>

<h4>Examples</h4>
<pre language="python">

</pre>
<hr>
<h2>plot_heat_power_density</h2>


        Function for plotting the heat flux profile



        <p><b>Returns</b>:

        None.

.
<h4>Parameters</h4>
<table cellpadding=10>
<tr><th>Name<th>Description
<tr><td><b><i>self</i></b><td>

</table>

<h4>Examples</h4>
<pre language="python">

</pre>
<hr>
<h2>set_coordinates</h2>


        Function for setting the local coordinates.



        <p><b>Parameters</b>:

        s_in : list or 1-by-n numpy array

            array of local coordinates.



        <p><b>Returns</b>:

        None.

.
<h4>Parameters</h4>
<table cellpadding=10>
<tr><th>Name<th>Description
<tr><td><b><i>self</i></b><td>s_in

</table>

<h4>Examples</h4>
<pre language="python">

</pre>
<hr>
<h2>set_edge_radius</h2>


        Function for setting the position of the high field side LCFS



        <p><b>Parameters</b>:

        radius_in : float

            position of the high field side LCFS.



        <p><b>Returns</b>:

        None.

.
<h4>Parameters</h4>
<table cellpadding=10>
<tr><th>Name<th>Description
<tr><td><b><i>self</i></b><td>radius_in

</table>

<h4>Examples</h4>
<pre language="python">

</pre>
<hr>
<h1 id="8ecd2f7bc5142547fcef6ceccbddd4a6">Class: InterfaceSurface</h1>
<p>

    Class for mapping power from an interface surface in the divertor onto the walls.



    :param Point2D point_a: A 2D point representing the start of the interface surface.

    :param Point2D point_b: A 2D point representing the end of the interface surface.

    :param ndarray power_profile: A an array of power values representing the power

      profile along the interface surface. These points are equally spaced and will be

      re-normalised to give an integral of one over the surface.</p>

<h2>__init__</h2>

<h4>Parameters</h4>
<table cellpadding=10>
<tr><th>Name<th>Description
<tr><td><b><i>self</i></b><td>point_a

</table>

<h4>Examples</h4>
<pre language="python">

</pre>
<hr>
<h2>_generate_sample_point</h2>

<h4>Parameters</h4>
<table cellpadding=10>
<tr><th>Name<th>Description
<tr><td><b><i>self</i></b><td>

</table>

<h4>Examples</h4>
<pre language="python">

</pre>
<hr>
<h2>_write_mesh_points_vtk</h2>

<h4>Parameters</h4>
<table cellpadding=10>
<tr><th>Name<th>Description
<tr><td><b><i>self</i></b><td>filename

</table>

<h4>Examples</h4>
<pre language="python">

</pre>
<hr>
<h2>_write_mesh_power_vtk</h2>

<h4>Parameters</h4>
<table cellpadding=10>
<tr><th>Name<th>Description
<tr><td><b><i>self</i></b><td>filename

</table>

<h4>Examples</h4>
<pre language="python">

</pre>
<hr>
<h2>get_lcfs_radius</h2>

<h4>Parameters</h4>
<table cellpadding=10>
<tr><th>Name<th>Description
<tr><td><b><i></i></b><td>

</table>

<h4>Examples</h4>
<pre language="python">

</pre>
<hr>
<h2>histogram_plot</h2>

<h4>Parameters</h4>
<table cellpadding=10>
<tr><th>Name<th>Description
<tr><td><b><i>self</i></b><td>

</table>

<h4>Examples</h4>
<pre language="python">

</pre>
<hr>
<h2>plot</h2>


        Plot the interface surface line across the equilibrium.



        :param EFITEquilibrium equilibrium: the equilibrium to plot.
<h4>Parameters</h4>
<table cellpadding=10>
<tr><th>Name<th>Description
<tr><td><b><i>EFITEquilibrium equilibrium</i></b><td>the equilibrium to plot.

</table>

<h4>Examples</h4>
<pre language="python">

</pre>
<hr>
<h2>psin</h2>

<h4>Parameters</h4>
<table cellpadding=10>
<tr><th>Name<th>Description
<tr><td><b><i>r</i></b><td>offset=0

</table>

<h4>Examples</h4>
<pre language="python">

</pre>
<hr>
<h1 id="4ab17cc669860a18b5efd64996345288">Class: ParticlePath()</h1>
<p>

    Class for determining the path of a charged particle in a magnetic field

    given a Fiesta equilibrium.



    Member functions:

        follow_path(p_0, max_time, max_points, rtol, break_at_limiter).</p>

<h2>__init__</h2>


        Class constructor, initialses variables and the interpolation functions

        to determine the B field, the gradient of magnitude of the B field, and

        the curl of the B field unit vector.
<h4>Parameters</h4>
<table cellpadding=10>
<tr><th>Name<th>Description
<tr><td><b><i>self</i></b><td>charge_mass_ratio

</table>

<h4>Examples</h4>
<pre language="python">

</pre>
<hr>
<h2>dlorentz_dt</h2>


            The function describing the ode to solve in order to track the

            location of patricle in a magnetic field (in accordance to the

            lorentz equation)



            <p><b>Input</b>: Time, time variable

                   Vec,  vector of the form [r, phi, z, v_r, v_phi, v_z]



            <p><b>Return</b>: DVecDt, the right-hand side of the ode to solve.
<h4>Parameters</h4>
<table cellpadding=10>
<tr><th>Name<th>Description
<tr><td><b><i>time</i></b><td>vec

</table>

<h4>Examples</h4>
<pre language="python">

</pre>
<hr>
<h2>dvec_dt</h2>


            The function describing the ode to solve in order to track the

            location of patricle in a magnetic field (in accordance with the

            guiding centre particle equation)



            <p><b>Input</b>: Time, time variable

                   Vec,  vector of the form [r, phi, z, v_para, moment]



            <p><b>Return</b>: DVecDt, the right-hand side of the ode to solve.
<h4>Parameters</h4>
<table cellpadding=10>
<tr><th>Name<th>Description
<tr><td><b><i>time</i></b><td>vec

</table>

<h4>Examples</h4>
<pre language="python">

</pre>
<hr>
<h2>event</h2>


            Function for determining whether the solution to the ODE passes a wall,

            which terminates the ODE solver



            <p><b>Input</b>: l_dist, np.array with the distance along the magnetic field-line

                   x_vec,  vector with the R, phi and Z initial positions



            <p><b>Return</b>: intersect_wall, returns 0 if any wall surface is intersected,

                                    otherwise returns a float (the event function of

                                    solve_ivp only looks for events = 0)

.
<h4>Parameters</h4>
<table cellpadding=10>
<tr><th>Name<th>Description
<tr><td><b><i>time</i></b><td>vec

</table>

<h4>Examples</h4>
<pre language="python">

</pre>
<hr>
<h2>event2</h2>


            Function for determining whether the solution to the ODE passes a wall,

            which terminates the ODE solver



            <p><b>Input</b>: l_dist, np.array with the distance along the magnetic field-line

                   x_vec,  vector with the R, phi and Z initial positions



            <p><b>Return</b>: intersect_wall, returns 0 if any wall surface is intersected,

                                    otherwise returns a float (the event function of

                                    solve_ivp only looks for events = 0).
<h4>Parameters</h4>
<table cellpadding=10>
<tr><th>Name<th>Description
<tr><td><b><i>time</i></b><td>vec

</table>

<h4>Examples</h4>
<pre language="python">

</pre>
<hr>
<h2>exact_b_phi</h2>


        calculate the exact Bphi value dependent on the radial position.
<h4>Parameters</h4>
<table cellpadding=10>
<tr><th>Name<th>Description
<tr><td><b><i>self</i></b><td>r_pos

</table>

<h4>Examples</h4>
<pre language="python">

</pre>
<hr>
<h1 id="2a84ce044b53bb3902f080782d296b09">Class: Vector2()</h1>
<p>Vectors in 2 dimensional euclidean space.</p>

<h2>Normalize</h2>
Normalize a copy of vector.
<h4>Parameters</h4>
<table cellpadding=10>
<tr><th>Name<th>Description
<tr><td><b><i>this</i></b><td>Vector to clone and normalize

</table>

<h4>Examples</h4>
<pre language="python">
assert((x * 3 + y * 4).Normalize() == x * 0.6 + y * 0.8)
</pre>
<hr>
<h2>R180</h2>
Rotate  a copy of a vector by 180 degrees.
<h4>Parameters</h4>
<table cellpadding=10>
<tr><th>Name<th>Description
<tr><td><b><i>this</i></b><td>Vector to clone and rotate

</table>

<h4>Examples</h4>
<pre language="python">
assert(x.R180() == -x)
</pre>
<hr>
<h2>R270</h2>
Rotate a vector by 270 degrees anticlockwise.
<h4>Parameters</h4>
<table cellpadding=10>
<tr><th>Name<th>Description
<tr><td><b><i>this</i></b><td>Vector to clone and rotate

</table>

<h4>Examples</h4>
<pre language="python">
assert(x.R270() == -y)
</pre>
<hr>
<h2>R90</h2>
Rotate a copy of a vector by 90 degrees anticlockwise.
<h4>Parameters</h4>
<table cellpadding=10>
<tr><th>Name<th>Description
<tr><td><b><i>this</i></b><td>Vector to clone and rotate

</table>

<h4>Examples</h4>
<pre language="python">
assert(x.R90() == y)
</pre>
<hr>
<h2>Swap</h2>
Swap the components of a copy of a vector.
<h4>Parameters</h4>
<table cellpadding=10>
<tr><th>Name<th>Description
<tr><td><b><i>this</i></b><td>Vector to clone and swap

</table>

<h4>Examples</h4>
<pre language="python">
z = x + y * 2
Z = z.Swap()
assert(z != Z)
</pre>
<hr>
<h2>__abs__</h2>
Length of a vector.
<h4>Parameters</h4>
<table cellpadding=10>
<tr><th>Name<th>Description
<tr><td><b><i>this</i></b><td>Vector whose length is required

</table>

<h4>Examples</h4>
<pre language="python">
assert(abs(x * 3 + y * 4) == 5)
</pre>
<hr>
<h2>__add__ **+ **</h2>
Add the second vector to a copy of the first vector.
<h4>Parameters</h4>
<table cellpadding=10>
<tr><th>Name<th>Description
<tr><td><b><i>this</i></b><td>Vector to be cloned and added to

<tr><td><b><i>other</i></b><td>Vector to add

</table>

<h4>Examples</h4>
<pre language="python">
assert(x + y == Vector2(1, 1))
</pre>
<hr>
<h2>__iadd__ **+=**</h2>
Add the second vector to the first vector.
<h4>Parameters</h4>
<table cellpadding=10>
<tr><th>Name<th>Description
<tr><td><b><i>this</i></b><td>Vector to be added to

<tr><td><b><i>other</i></b><td>Vector to add

</table>

<h4>Examples</h4>
<pre language="python">
X = x
X = X + x
assert(x == Vector2(1, 0) and X == Vector2(2, 0))
</pre>
<hr>
<h2>__imul__ ***=**</h2>
Multiply a vector by a scalar.
<h4>Parameters</h4>
<table cellpadding=10>
<tr><th>Name<th>Description
<tr><td><b><i>this</i></b><td>Vector to be multiplied

<tr><td><b><i>scalar</i></b><td>Scalar to multiply by

</table>

<h4>Examples</h4>
<pre language="python">
X = x
X = X * 2
assert(x == Vector2(1, 0) and X == Vector2(2, 0))
</pre>
<hr>
<h2>__init__</h2>
Create a vector.
<h4>Parameters</h4>
<table cellpadding=10>
<tr><th>Name<th>Description
<tr><td><b><i>this</i></b><td>New vector to be initialized

<tr><td><b><i>x = 0</i></b><td>X coordinate

<tr><td><b><i>y = 0</i></b><td>Y coordinate

</table>

<h4>Examples</h4>
<pre language="python">
assert(x + y == Vector2(1, 1))
</pre>
<hr>
<h2>__isub__ **-=**</h2>
Subtract the second vector from the first vector.
<h4>Parameters</h4>
<table cellpadding=10>
<tr><th>Name<th>Description
<tr><td><b><i>this</i></b><td>Vector to be subtracted from

<tr><td><b><i>other</i></b><td>Vector to subtract

</table>

<h4>Examples</h4>
<pre language="python">
X = x
X = X - x
assert(x == Vector2(1, 0) and X == Vector2(0, 0))
</pre>
<hr>
<h2>__itruediv__ **/=**</h2>
Divide a vector by a scalar.
<h4>Parameters</h4>
<table cellpadding=10>
<tr><th>Name<th>Description
<tr><td><b><i>this</i></b><td>Vector to be divided

<tr><td><b><i>scalar</i></b><td>Scalar to divide by

</table>

<h4>Examples</h4>
<pre language="python">
X = x * 4
X = X / 2
assert(x == Vector2(1, 0) and X == Vector2(2, 0))
</pre>
<hr>
<h2>__mul__ *** **</h2>
Multiply a copy of vector by a scalar.
<h4>Parameters</h4>
<table cellpadding=10>
<tr><th>Name<th>Description
<tr><td><b><i>this</i></b><td>Vector to be cloned and multiplied

<tr><td><b><i>scalar</i></b><td>Scalar to multiply by

</table>

<h4>Examples</h4>
<pre language="python">
assert(x * 2 + y * 3 == Vector2(2, 3))
</pre>
<hr>
<h2>__neg__ **- **</h2>
Rotate a copy of a vector by 180 degrees.
<h4>Parameters</h4>
<table cellpadding=10>
<tr><th>Name<th>Description
<tr><td><b><i>this</i></b><td>Vector to be cloned and negated

</table>

<h4>Examples</h4>
<pre language="python">
assert(- y == Vector2(0, -1))
</pre>
<hr>
<h2>__repr__</h2>
String representation of a vector.
<h4>Parameters</h4>
<table cellpadding=10>
<tr><th>Name<th>Description
<tr><td><b><i>this</i></b><td>Vector to be represented

</table>

<h4>Examples</h4>
<pre language="python">
assert(repr(x + y * 2) == 'Vector2(1, 2)')
</pre>
<hr>
<h2>__sub__ **- **</h2>
Subtract the second vector from a copy of the first vector.
<h4>Parameters</h4>
<table cellpadding=10>
<tr><th>Name<th>Description
<tr><td><b><i>this</i></b><td>Vector to be cloned and subtracted from

<tr><td><b><i>other</i></b><td>Vector to subtract

</table>

<h4>Examples</h4>
<pre language="python">
assert(x - y == Vector2(1, -1))
</pre>
<hr>
<h2>__truediv__ **/ **</h2>
Divide a copy of a vector by a scalar.
<h4>Parameters</h4>
<table cellpadding=10>
<tr><th>Name<th>Description
<tr><td><b><i>this</i></b><td>Vector to be cloned and divided

<tr><td><b><i>scalar</i></b><td>Scalar to divide by

</table>

<h4>Examples</h4>
<pre language="python">
assert((x * 4 + y * 2) / 2 == x * 2 + y)
</pre>
<hr>
<h2>angle</h2>
Angle in radians anticlockwise that the first vector must be rotated to point along the second vector normalized to the range: -pi to +pi.
<h4>Parameters</h4>
<table cellpadding=10>
<tr><th>Name<th>Description
<tr><td><b><i>o</i></b><td>Left hand vector we are looking along

<tr><td><b><i>p</i></b><td>Right hand vector we want anticlockwise angle to

</table>

<h4>Examples</h4>
<pre language="python">
for i in range(-179, +179): # Anticlockwise angle from x
assert(Vector2.close
(x.angle(x * math.cos(dr(i)) + y * math.sin(dr(i))), dr(i)))
</pre>
<hr>
<h2>area</h2>
Signed area of the parallelogram defined by the two vectors. The area is negative if the second vector appears to the right of the first if they are both placed at the origin and the observer stands against the z-axis in a left handed coordinate system.
<h4>Parameters</h4>
<table cellpadding=10>
<tr><th>Name<th>Description
<tr><td><b><i>this</i></b><td>Left hand vector

<tr><td><b><i>other</i></b><td>Right hand vector

</table>

<h4>Examples</h4>
<pre language="python">
assert((x + y).area(-x + y) == 2)
</pre>
<hr>
<h2>clone</h2>
Clone a vector to allow it to be modified by other operations.
<h4>Parameters</h4>
<table cellpadding=10>
<tr><th>Name<th>Description
<tr><td><b><i>this</i></b><td>Vector to be cloned

</table>

<h4>Examples</h4>
<pre language="python">
z = x + y * 2
Z = z.clone()
assert(z == Z)
</pre>
<hr>
<h2>close</h2>
Whether two numbers are close.
<h4>Parameters</h4>
<table cellpadding=10>
<tr><th>Name<th>Description
<tr><td><b><i>a</i></b><td>First number to compare

<tr><td><b><i>b</i></b><td>Second number to compare

</table>

<h4>Examples</h4>
<pre language="python">
assert(Vector2.close(0, Vector2.closeRadius() / 2))
</pre>
<hr>
<h2>cos</h2>
cos(angle between two vectors).
<h4>Parameters</h4>
<table cellpadding=10>
<tr><th>Name<th>Description
<tr><td><b><i>this</i></b><td>First vector

<tr><td><b><i>other</i></b><td>Second vector

</table>

<h4>Examples</h4>
<pre language="python">
assert(Vector2.close((x + y).cos(y), 1 / r2))
assert(Vector2.close( x.cos(x + yr3), 0.5))
</pre>
<hr>
<h2>distance</h2>
Distance between the points identified by two vectors when placed on the same point.
<h4>Parameters</h4>
<table cellpadding=10>
<tr><th>Name<th>Description
<tr><td><b><i>this</i></b><td>Vector to start point

<tr><td><b><i>other</i></b><td>Vector to end point

</table>

<h4>Examples</h4>
<pre language="python">
assert((x * 3 + y * 4).distance (-(x * 3 + y * 4)) == 10)
</pre>
<hr>
<h2>distance2</h2>
Distance squared between the points identified

       by two vectors when placed on the same point.
<h4>Parameters</h4>
<table cellpadding=10>
<tr><th>Name<th>Description
<tr><td><b><i>this</i></b><td>Vector to start point

<tr><td><b><i>other</i></b><td>Vector to end point

</table>

<h4>Examples</h4>
<pre language="python">
assert((x * 3 + y * 4).distance2(-(x * 3 + y * 4)) == 100)
</pre>
<hr>
<h2>dot</h2>
Dot product of two vectors.
<h4>Parameters</h4>
<table cellpadding=10>
<tr><th>Name<th>Description
<tr><td><b><i>this</i></b><td>First vector

<tr><td><b><i>other</i></b><td>Second vector

</table>

<h4>Examples</h4>
<pre language="python">
assert((x * 2 + y).dot(x + y * 3) == 5)
</pre>
<hr>
<h2>length</h2>
Length of a vector.
<h4>Parameters</h4>
<table cellpadding=10>
<tr><th>Name<th>Description
<tr><td><b><i>this</i></b><td>Vector whose length is required

</table>

<h4>Examples</h4>
<pre language="python">
assert((x * 3 + y * 4).length() == 5)
</pre>
<hr>
<h2>length2</h2>
Length squared of a vector.
<h4>Parameters</h4>
<table cellpadding=10>
<tr><th>Name<th>Description
<tr><td><b><i>this</i></b><td>Vectors whose length squared is required

</table>

<h4>Examples</h4>
<pre language="python">
assert((x * 3 + y * 4).length2() == 25)
</pre>
<hr>
<h2>normalize</h2>
Normalize a vector.
<h4>Parameters</h4>
<table cellpadding=10>
<tr><th>Name<th>Description
<tr><td><b><i>this</i></b><td>Vector to normalize

</table>

<h4>Examples</h4>
<pre language="python">
assert((x * 3 + y * 4).clone().normalize() == x * 0.6 + y * 0.8)
</pre>
<hr>
<h2>r180</h2>
Rotate a vector by 180 degrees.
<h4>Parameters</h4>
<table cellpadding=10>
<tr><th>Name<th>Description
<tr><td><b><i>this</i></b><td>Vector to rotate

</table>

<h4>Examples</h4>
<pre language="python">
assert(x.clone().r180() == -x)
</pre>
<hr>
<h2>r270</h2>
Rotate a copy of a vector by 270 degrees anticlockwise.
<h4>Parameters</h4>
<table cellpadding=10>
<tr><th>Name<th>Description
<tr><td><b><i>this</i></b><td>Vector to rotate

</table>

<h4>Examples</h4>
<pre language="python">
assert(x.clone().r270() == -y)
</pre>
<hr>
<h2>r90</h2>
Rotate a vector by 90 degrees anticlockwise.
<h4>Parameters</h4>
<table cellpadding=10>
<tr><th>Name<th>Description
<tr><td><b><i>this</i></b><td>Vector to rotate

</table>

<h4>Examples</h4>
<pre language="python">
assert(x.clone().r90() == y)
</pre>
<hr>
<h2>sin</h2>
sin(angle between two vectors).
<h4>Parameters</h4>
<table cellpadding=10>
<tr><th>Name<th>Description
<tr><td><b><i>this</i></b><td>Left hand vector

<tr><td><b><i>other</i></b><td>Right hand vector

</table>

<h4>Examples</h4>
<pre language="python">
assert(Vector2.close((x + y).sin(y), 1 / r2))
assert(Vector2.close( x.sin(x + yr3), r3 / 2))
</pre>
<hr>
<h2>smallestAngleToNormalPlane</h2>
The smallest angle between the second vector and a plane normal to the first vector.
<h4>Parameters</h4>
<table cellpadding=10>
<tr><th>Name<th>Description
<tr><td><b><i>a</i></b><td>Vector normal to plane

<tr><td><b><i>b</i></b><td>Vector at angle to plane

</table>

<h4>Examples</h4>
<pre language="python">
assert(Vector2.close(dr( 0), y.smallestAngleToNormalPlane( x))) # First vector is y, second vector is 0 degrees anti-clockwise from x axis
assert(Vector2.close(dr(+45), y.smallestAngleToNormalPlane( x + y))) # +45
assert(Vector2.close(dr(+90), y.smallestAngleToNormalPlane( y))) # +90
assert(Vector2.close(dr(+45), y.smallestAngleToNormalPlane(-x + -y))) # +135
assert(Vector2.close(dr( 0), y.smallestAngleToNormalPlane(-x))) # +180
assert(Vector2.close(dr(+45), y.smallestAngleToNormalPlane(-x + -y))) # +225
assert(Vector2.close(dr(+90), y.smallestAngleToNormalPlane( -y))) # +270
assert(Vector2.close(dr(+45), y.smallestAngleToNormalPlane(-x + -y))) # +315
assert(Vector2.close(dr( 0), y.smallestAngleToNormalPlane( x))) # +360
</pre>
<hr>
<h2>swap</h2>
Swap the components of a vector.
<h4>Parameters</h4>
<table cellpadding=10>
<tr><th>Name<th>Description
<tr><td><b><i>this</i></b><td>Vector to swap

</table>

<h4>Examples</h4>
<pre language="python">
assert((x + y * 2).swap() == y + x * 2)
</pre>
<hr>
<h1 id="81c15fb95b1c3422ceb5b2b85c7f9668">Class: add_resource.py</h1>
<p></p>

<h2>main</h2>

<h4>Parameters</h4>
<table cellpadding=10>
<tr><th>Name<th>Description
<tr><td><b><i></i></b><td>

</table>

<h4>Examples</h4>
<pre language="python">

</pre>
<hr>
<h1 id="afea94043c155f8b6a9c8265b54802f8">Class: calculate_angle.py</h1>
<p></p>

<h2>_angle</h2>


        Calculate the angle between two vectors defined as:



            cos(angle) = (v_1 dot v_2)/(||v_1|| * ||v_2||),



        where ||v|| is the length of the vector.



        <p><b>Parameters</b>:

        v_1 : 2-x-1 np.array

            Array with x and y of the first vector

        v_2 : 2-x-1 np.array

            Array with x and y of the second vector



        <p><b>Returns</b>:

        angle : float

            A float with the angle between the two vectors

.
<h4>Parameters</h4>
<table cellpadding=10>
<tr><th>Name<th>Description
<tr><td><b><i>v_1</i></b><td>v_2

</table>

<h4>Examples</h4>
<pre language="python">

</pre>
<hr>
<h2>calculate_angle</h2>


    Function for calculating the angle with respect to the normal of v_1

    between two vectors.



    <p><b>Parameters</b>:

    v_1_vectors : 2-x-1 np.array

        2D vector of the form (x, y)

    v_2 : 2-x-1 np.array

        2D vector of the form (x, y)



    <p><b>Returns</b>:

    angle : float

        The angle between the two vectors

.
<h4>Parameters</h4>
<table cellpadding=10>
<tr><th>Name<th>Description
<tr><td><b><i>v_1</i></b><td>v_2

</table>

<h4>Examples</h4>
<pre language="python">

</pre>
<hr>
<h1 id="2bc2db0b2b880fd2ffc72ac6bdb18365">Class: equdsk.py</h1>
<p></p>

<h2>read_eqdsk</h2>

<h4>Parameters</h4>
<table cellpadding=10>
<tr><th>Name<th>Description
<tr><td><b><i></i></b><td>

</table>

<h4>Examples</h4>
<pre language="python">

</pre>
<hr>
<h1 id="f7f0f2db6aa5ccab3d80ff8536dbf6f7">Class: field_line_projection.py</h1>
<p></p>

<h2>_flux_expansion</h2>


    Function for evaluating the flux expansion, defined as:



        f_x = x_omp*b_pol(x_omp, y_omp)/(x_div*b_pol(x_div, y_div)),



    where x_omp is the radial coordinate at the OMP, y_omp is the vertical

    coordinate at the OMP (usually 0) and x_div and y_div are the corresponding

    coordinates at the divertor.



    <p><b>Parameters</b>:

    b_pol : cherab.core.math.Interpolate2DCubic

        Function for interpolating the poloidal magnetic field from a fiesta

        equilibrium, given an (r, z)-point

    points_omp : 2-x-1 list

        List containing a point at the outboard mid-plane

    points_div : 2-x-1 list

        List containing the point at the corresponding flux surface at the

        divertor.



    <p><b>Returns</b>:

    f_x : float

        The flux expansion at the divertor point specified

.
<h4>Parameters</h4>
<table cellpadding=10>
<tr><th>Name<th>Description
<tr><td><b><i>b_pol</i></b><td>points_omp

</table>

<h4>Examples</h4>
<pre language="python">

</pre>
<hr>
<h2>project_field_lines</h2>


    Function mapping the field-lines from the specified coordinates at the

    OMP to the specified coordinates at a given surface. Currently the surface

    is assumed to be represented by a 1D polynomial function, y = ax + b.



    <p><b>Parameters</b>:

    x_axis_omp : n-x-1 np.array

        Numpy array with the radial coordinates we wish to map at the OMP

    fiesta : Fiesta

        A Fiesta object with the 2D equilibrium we wish to map

    divertor_coords : 2-x-2 np.array

        A 2-x-2 numpy array containg the corner points of the divertor in the

        2D projection



    <p><b>Returns</b>:

    divertor_map : dictionary

        A dictionary containing:

            "R_div" : an n-x-1 array

                with the R-coordinates at the divertor tile

                corresponding to the same psi_n as at the OMP

            "Z_div" : an n-x-1 array

                with the Z-coordinates at the divertor tile

                corresponding to the same psi_n as at the OMP

            "Angles" : an n-x-1 array

                with the angles between the field lines and the divertor tile

                corresponding to the same psi_n as at the OMP

            "Flux_expansion" : an n-x-1 array

                with the flux expasion at the divertor tile

                corresponding to the same psi_n as at the OMP

.
<h4>Parameters</h4>
<table cellpadding=10>
<tr><th>Name<th>Description
<tr><td><b><i>x_axis_omp</i></b><td>surface_coords

</table>

<h4>Examples</h4>
<pre language="python">

</pre>
<hr>
<h1 id="431e5d48b3ab6a117ee828bfec778106">Class: getOption.py</h1>
<p></p>

<h2>getOption</h2>
Get the named long option from the command line without modifying the command line arguments array.
<h4>Parameters</h4>
<table cellpadding=10>
<tr><th>Name<th>Description
<tr><td><b><i>option</i></b><td>

</table>

<h4>Examples</h4>
<pre language="python">

</pre>
<hr>
<h1 id="9386ee3547cc05f3bdb7852a4e7fbdfa">Class: intersections.py</h1>
<p></p>

<h2>_get_rectangle_intersections</h2>


    Function dividing two functions of [x,y] coordinates into rectangles

    corresponding to the number of elements in each function and evaluating

    the indices where the rectangles intersect.



    <p><b>Input</b>: func1, a numpy array with the two numpy arrays corresponding to x and y

                  for the first function

           func2, a numpy array with the two numpy arrays corresponding to x and y

                  for the second function



    <p><b>Return</b>: (i, j), a tuple where

                    i is a numpy array with the indices for the

                    intersections in the first function and

                    j is a numpy array with the indices for the

                    intersections in the second function.
<h4>Parameters</h4>
<table cellpadding=10>
<tr><th>Name<th>Description
<tr><td><b><i>func1</i></b><td>func2

</table>

<h4>Examples</h4>
<pre language="python">

</pre>
<hr>
<h2>intersection</h2>


    Function for calculated the intersection between two curves.

    Computes the (x,y) locations where two curves intersect.



    The theory is;

    Given two line segments, L1 and L2,



    with L1 endpoints:  (x1(1),y1(1)) and (x1(2),y1(2))

    and  L2 endpoints:  (x2(1),y2(1)) and (x2(2),y2(2))



    we can write four equations with four unknowns and then solve them.  The

    four unknowns are t1, t2, x0 and y0, where (x0,y0) is the intersection of

    L1 and L2, t1 is the distance from the starting point of L1 to the

    intersection relative to the length of L1 and t2 is the distance from the

    starting point of L2 to the intersection relative to the length of L2.

    So, the four equations are



        (x1(2) - x1(1))*t1 = x0 - x1(1)

        (x2(2) - x2(1))*t2 = x0 - x2(1)

        (y1(2) - y1(1))*t1 = y0 - y1(1)

        (y2(2) - y2(1))*t2 = y0 - y2(1)



    Rearranging and writing in matrix form gives



        [x1(2)-x1(1)       0       -1   0;      [t1;      [-x1(1);

              0       x2(2)-x2(1)  -1   0;   *   t2;   =   -x2(1);

         y1(2)-y1(1)       0        0  -1;       x0;       -y1(1);

              0       y2(2)-y2(1)   0  -1]       y0]       -y2(1)]



    Let's call that A*T = B.  We can solve for T with T = A/B.



    Once we have our solution we just have to look at t1 and t2 to determine

    whether L1 and L2 intersect.  If 0 <= t1 < 1 and 0 <= t2 < 1 then the two

    line segments cross and we can include (x0,y0) in the output.



    To avoid having to do this for every line segment, it is checked if the line

    segments can possibly intersect by dividing line segments into rectangles

    and testing for an overlap between the triangles.



    <p><b>Input</b>: func1, a numpy array with the two numpy arrays corresponding to x and y

                  for the first function

           func2, a numpy array with the two numpy arrays corresponding to x and y

                  for the second function



    <p><b>Return</b>: i,    a numpy array of floats with the sum of the indices and distances

                  [0; 1[ to the intersections of func1

            j,    a numpy array of floats with the sum of the indices and distances

                  [0; 1[ to the intersections of func2

            x0,   a numpy array with the x positions of the intersections

            y0,   a numpy array with the y positions of the intersections.
<h4>Parameters</h4>
<table cellpadding=10>
<tr><th>Name<th>Description
<tr><td><b><i>func1</i></b><td>func2

</table>

<h4>Examples</h4>
<pre language="python">

</pre>
<hr>
<h1 id="b53b675887a97eff5d9fe276b5c074bb">Class: load_pickle.py</h1>
<p></p>

<h2>load_pickle</h2>


    Function for loading a .pickle file



    <p><b>Input</b>: file_name,    a string with the .pickle file to load, e.g. 'fiesta/eq_0002'



    <p><b>Return</b>: pickle_dict, the dictionary stored in the .pickle file.
<h4>Parameters</h4>
<table cellpadding=10>
<tr><th>Name<th>Description
<tr><td><b><i>file_name</i></b><td>

</table>

<h4>Examples</h4>
<pre language="python">

</pre>
<hr>
<h1 id="e02daeaf1215d7e93399246760393699">Class: machine_configuration.py</h1>
<p></p>

<h2>load_wall_configuration</h2>

<h4>Parameters</h4>
<table cellpadding=10>
<tr><th>Name<th>Description
<tr><td><b><i>config_file</i></b><td>parent

</table>

<h4>Examples</h4>
<pre language="python">

</pre>
<hr>
<h1 id="225106ddf5f55966717a2e58dbd07ea3">Class: map_field_lines.py</h1>
<p></p>

<h2>map_field_lines</h2>


    Function for mapping each field-line to the intersection with the vessel walls



    <p><b>Input</b>: x_vec_at_omp,  a numpy array with the radial points at the OMP where

                          we want the mapping to start

           file_path,     a string with the path to the FIESTA equilibrium .mat file

           configuration, a string with either 'limited' or 'diverted', where 'diverted'

                          is the defualt configuration



    <p><b>Return</b>: field_line_dict, a python dictionary with the radial position from the OMP in m

                             as the key and the field-line dictionary with the R, phi

                             and Z components along the field line as well as the length, l,

                             from the LCFS to the current point along the field-line.
<h4>Parameters</h4>
<table cellpadding=10>
<tr><th>Name<th>Description
<tr><td><b><i>x_vec_at_omp</i></b><td>file_path

</table>

<h4>Examples</h4>
<pre language="python">

</pre>
<hr>
<h1 id="af5a643bca05127c7079d216ab5fd94c">Class: midplane_power.py</h1>
<p></p>

<h2>run_midplane_power</h2>


    Function for running the specified mid-plane model



    <p><b>Input</b>: midplane_model, a string with the type of model to use for the heat-flux

                           evaluation

           plasma,         the plasma settings defined in the .json input file



    <p><b>Return</b>: footprint,     an object of the midplane_model class specified.
<h4>Parameters</h4>
<table cellpadding=10>
<tr><th>Name<th>Description
<tr><td><b><i>midplane_model</i></b><td>plasma

</table>

<h4>Examples</h4>
<pre language="python">

</pre>
<hr>
<h1 id="56845d841c468ee501c30899177b17c6">Class: midplane_to_divertor_power.py</h1>
<p></p>

<h2>run_midplane_to_div_power</h2>


    .
<h4>Parameters</h4>
<table cellpadding=10>
<tr><th>Name<th>Description
<tr><td><b><i>midplane_model</i></b><td>plasma

</table>

<h4>Examples</h4>
<pre language="python">

</pre>
<hr>
<h1 id="9fb2a57520e98a6c2e8897b9609268d9">Class: parser.py</h1>
<p></p>

<h2>is_valid_file</h2>


    Function for checking it a file exists



    <p><b>Input</b>:  parser, an argparse.ArgumentParser object for parsing the terminal input

            arg,    a terminal input with the file to load



    <p><b>Return</b>: open(arg, 'r'), an open file handle if file exists, else -1.
<h4>Parameters</h4>
<table cellpadding=10>
<tr><th>Name<th>Description
<tr><td><b><i>parser</i></b><td>arg

</table>

<h4>Examples</h4>
<pre language="python">

</pre>
<hr>
<h1 id="c49c1165dd9fa01598d89310eece1054">Class: physics_constants.py</h1>
<p></p>

<h2>get_physics_constants</h2>


    Function for populating the physics constant dictionary



    <p><b>Input</b>: none



    <p><b>Return</b>: physics_constants, a dictionary with a series of physics constants

            in SI units.
<h4>Parameters</h4>
<table cellpadding=10>
<tr><th>Name<th>Description
<tr><td><b><i></i></b><td>

</table>

<h4>Examples</h4>
<pre language="python">

</pre>
<hr>
<h1 id="cc87324b013b0c555816970d52f4e4d4">Class: project_heat_flux.py</h1>
<p></p>

<h2>project_heat_flux</h2>


    Function for mapping the heat flux from the OMP to the divertor. The heat

    flux at a different position is given by:



        q_parallel_surf = R_omp/R_surf*q_parallel_omp/(f_x/cos(alpha)),



    where R_omp, is the radial coordinate at the OMP, R_surf is the radial coordinates

    at the surface, q_parallel_omp is the parallel heat flux at the OMP, alpha

    is the incidence angle of the field-lines with respect to the normal of the

    surface, and f_x is the flux expansion:

        

        f_x = R_omp*B_pol(R_omp, Z_omp)/(R_surf*B_pol(R_surf, Z_surf)),



    where Z_omp is the vertical position of the OMP (usually 0), Z_surf is

    the vertical position of the surface, and B_pol is the poloidal magnetic

    field and the given coordinates.



    <p><b>Parameters</b>:

    x_pos_omp : n-by-1 np.array

        Radial coordinates at the OMP

    heat_flux_profile : n-by-1 np.array

        Parallel heat flux at the given coordinates

    map_dict : dictionary

        Python dictionary with:

            keys: float

                x_pos_omp[i], each position at the omp has a corresponding mapped position

            values: dictionary

                dictionary with keys:

                    "R_pos" : float, radial position of the surface

                    "Z_pos" : float, vertical position of the surface

                    "f_x"   : float, flux expansion at the given R, Z position

                    "alpha" : float, incidence angle with respect to the normal

                                     of the surface

                



    <p><b>Returns</b>:

    q_surf : n-by-1 np.array

        parallel heat flux at the surface position

.
<h4>Parameters</h4>
<table cellpadding=10>
<tr><th>Name<th>Description
<tr><td><b><i>x_pos_omp</i></b><td>heat_flux_profile

</table>

<h4>Examples</h4>
<pre language="python">

</pre>
<hr>
<h1 id="537ca69b729227582c80ebad18406ae1">Class: psi_map_projection.py</h1>
<p></p>

<h2>_angle</h2>


        Calculate the angle between two vectors defined as:



            cos(angle) = (v_1 dot v_2)/(||v_1|| * ||v_2||),



        where ||v|| is the length of the vector.



        <p><b>Parameters</b>:

        v1 : 2-x-1 np.array

            Array with x and y of the first vector

        v2 : 2-x-1 np.array

            Array with x and y of the first vector



        <p><b>Returns</b>:

        angle : float

            A float with the angle between the two vectors

.
<h4>Parameters</h4>
<table cellpadding=10>
<tr><th>Name<th>Description
<tr><td><b><i>v_1</i></b><td>v_2

</table>

<h4>Examples</h4>
<pre language="python">

</pre>
<hr>
<h2>_calculate_angles</h2>


    Function for calculating the angles between the normalised psi calculated

    at the divertor and the divertor itself.

    The angle specified is calculated as the magnitude of the angle between the

    normal of the divertor and the normalised psi contour.



    <p><b>Parameters</b>:

    v_1_vectors : n-x-1 list of 2-x-1 np.arrays

        List of vectors of the field-lines close to the divertor interface

    v_2 : 2-x-1 np.array

        Vector for the divertor coordinates



    <p><b>Returns</b>:

    angles : n-x-1 np.array

        Array of angles between the divertor and the normalised psi calculated

        at the divertor.

.
<h4>Parameters</h4>
<table cellpadding=10>
<tr><th>Name<th>Description
<tr><td><b><i>v_1_vectors</i></b><td>v_2

</table>

<h4>Examples</h4>
<pre language="python">

</pre>
<hr>
<h2>_flux_expansion</h2>


    Function for evaluating the flux expansion, defined as:



        f_x = x_omp*b_pol(x_omp, y_omp)/(x_div*b_pol(x_div, y_div)),



    where x_omp is the radial coordinate at the OMP, y_omp is the vertical

    coordinate at the OMP (usually 0) and x_div and y_div are the corresponding

    coordinates at the divertor.



    <p><b>Parameters</b>:

    b_pol : cherab.core.math.Interpolate2DCubic

        Function for interpolating the poloidal magnetic field from a fiesta

        equilibrium, given an (r, z)-point

    points_omp : 2-x-1 list

        List containing a point at the outboard mid-plane

    points_div : 2-x-1 list

        List containing the point at the corresponding flux surface at the

        divertor.



    <p><b>Returns</b>:

    f_x : float

        The flux expansion at the divertor point specified

.
<h4>Parameters</h4>
<table cellpadding=10>
<tr><th>Name<th>Description
<tr><td><b><i>b_pol</i></b><td>points_omp

</table>

<h4>Examples</h4>
<pre language="python">

</pre>
<hr>
<h2>map_psi_omp_to_divertor</h2>


    Function mapping the normalised psi from the specified coordinates at the

    OMP to the specified coordinates at the divertor. Currently the divertor is

    assumed to be represented by a 1D polynomial function, y = ax + b.



    <p><b>Parameters</b>:

    x_axis_omp : n-x-1 np.array

        Numpy array with the radial coordinates we wish to map at the OMP

    fiesta : Fiesta

        A Fiesta object with the 2D equilibrium we wish to map

    divertor_coords : 2-x-2 np.array

        A 2-x-2 numpy array containg the corner points of the divertor in the

        2D projection



    <p><b>Returns</b>:

    divertor_map : dictionary

        A dictionary containing:

            "R_div" : an n-x-1 array

                with the R-coordinates at the divertor tile

                corresponding to the same psi_n as at the OMP

            "Z_div" : an n-x-1 array

                with the Z-coordinates at the divertor tile

                corresponding to the same psi_n as at the OMP

            "Angles" : an n-x-1 array

                with the angles between the field lines and the divertor tile

                corresponding to the same psi_n as at the OMP

            "Flux_expansion" : an n-x-1 array

                with the flux expasion at the divertor tile

                corresponding to the same psi_n as at the OMP

.
<h4>Parameters</h4>
<table cellpadding=10>
<tr><th>Name<th>Description
<tr><td><b><i>x_axis_omp</i></b><td>divertor_coords

</table>

<h4>Examples</h4>
<pre language="python">

</pre>
<hr>
<h1 id="4763e95919953c2a5b089bdd3fafb4d5">Class: remove_resource.py</h1>
<p></p>

<h2>main</h2>

<h4>Parameters</h4>
<table cellpadding=10>
<tr><th>Name<th>Description
<tr><td><b><i></i></b><td>

</table>

<h4>Examples</h4>
<pre language="python">

</pre>
<hr>
<h1 id="7db598b19ebf90a0f9549002f7f2d962">Class: resource_manager.py</h1>
<p></p>

<h2>_test_allowed_characters</h2>

<h4>Parameters</h4>
<table cellpadding=10>
<tr><th>Name<th>Description
<tr><td><b><i>string</i></b><td>

</table>

<h4>Examples</h4>
<pre language="python">

</pre>
<hr>
<h2>add_resource</h2>


    Add a resource file to the VITA hidden resource directory.



    :param str machine: A string identifying the fusion machine (e.g. "ST40").

    :param str type: The type of resource file. Must be one of ["equilibrium", "mesh"]

    :param str id: A unique string identifier for the resource.

    :param str path: The file path to the selected resource.

    :param bool symlink: Flag to specify whether a symlink should be created instead of

      copying the file. Useful for large simulation files.
<h4>Parameters</h4>
<table cellpadding=10>
<tr><th>Name<th>Description
<tr><td><b><i>str machine</i></b><td>A string identifying the fusion machine (e.g. "ST40").

<tr><td><b><i>str type</i></b><td>The type of resource file. Must be one of ["equilibrium", "mesh"]

<tr><td><b><i>str id</i></b><td>A unique string identifier for the resource.

<tr><td><b><i>str path</i></b><td>The file path to the selected resource.

<tr><td><b><i>bool symlink</i></b><td>Flag to specify whether a symlink should be created instead of

</table>

<h4>Examples</h4>
<pre language="python">

</pre>
<hr>
<h2>get_resource</h2>


    Get the path to a VITA resource.



    :param str machine: A string identifying the fusion machine (e.g. "ST40").

    :param str type: The type of resource file. Must be one of ["equilibrium", "mesh"]

    :param str id: A unique string identifier for the resource.

    :<p><b>Return</b>: The path to the selected resource.
<h4>Parameters</h4>
<table cellpadding=10>
<tr><th>Name<th>Description
<tr><td><b><i>str machine</i></b><td>A string identifying the fusion machine (e.g. "ST40").

<tr><td><b><i>str type</i></b><td>The type of resource file. Must be one of ["equilibrium", "mesh"]

<tr><td><b><i>str id</i></b><td>A unique string identifier for the resource.

</table>

<h4>Examples</h4>
<pre language="python">

</pre>
<hr>
<h2>list_resources</h2>

<h4>Parameters</h4>
<table cellpadding=10>
<tr><th>Name<th>Description
<tr><td><b><i></i></b><td>

</table>

<h4>Examples</h4>
<pre language="python">

</pre>
<hr>
<h2>remove_resource</h2>


    Remove a resource from the VITA resource file catalog.



    :param str machine: A string identifying the fusion machine (e.g. "ST40").

    :param str type: The type of resource file. Must be one of ["equilibrium", "mesh"]

    :param str id: A unique string identifier for the resource.

    :param bool prompt_user: Ask the user to confirm before deleting the file.
<h4>Parameters</h4>
<table cellpadding=10>
<tr><th>Name<th>Description
<tr><td><b><i>str machine</i></b><td>A string identifying the fusion machine (e.g. "ST40").

<tr><td><b><i>str type</i></b><td>The type of resource file. Must be one of ["equilibrium", "mesh"]

<tr><td><b><i>str id</i></b><td>A unique string identifier for the resource.

<tr><td><b><i>bool prompt_user</i></b><td>Ask the user to confirm before deleting the file.

</table>

<h4>Examples</h4>
<pre language="python">

</pre>
<hr>
<h2>update_resource</h2>


    Over write an existing resource file in the VITA hidden resource directory.



    :param str machine: A string identifying the fusion machine (e.g. "ST40").

    :param str type: The type of resource file. Must be one of ["equilibrium", "mesh"]

    :param str id: A unique string identifier for the resource.

    :param str path: The file path to the selected resource.

    :param bool symlink: Flag to specify whether a symlink should be created instead of

      copying the file. Useful for large simulation files.
<h4>Parameters</h4>
<table cellpadding=10>
<tr><th>Name<th>Description
<tr><td><b><i>str machine</i></b><td>A string identifying the fusion machine (e.g. "ST40").

<tr><td><b><i>str type</i></b><td>The type of resource file. Must be one of ["equilibrium", "mesh"]

<tr><td><b><i>str id</i></b><td>A unique string identifier for the resource.

<tr><td><b><i>str path</i></b><td>The file path to the selected resource.

<tr><td><b><i>bool symlink</i></b><td>Flag to specify whether a symlink should be created instead of

</table>

<h4>Examples</h4>
<pre language="python">

</pre>
<hr>
<h1 id="6612d973c9ef53d7152af2c08956fc7a">Class: update_resource.py</h1>
<p></p>

<h2>main</h2>

<h4>Parameters</h4>
<table cellpadding=10>
<tr><th>Name<th>Description
<tr><td><b><i></i></b><td>

</table>

<h4>Examples</h4>
<pre language="python">

</pre>
<hr>
<h1 id="9ebf128f7c57b374928a8d28e3bc06f4">Class: vitarunner.py</h1>
<p></p>

<h2>main</h2>

<h4>Parameters</h4>
<table cellpadding=10>
<tr><th>Name<th>Description
<tr><td><b><i></i></b><td>

</table>

<h4>Examples</h4>
<pre language="python">

</pre>
<hr>
<h1>Possible improvements to documentation</h1>
<h2>/home/phil/vita/core/vita/controller/midplane_power.py</h2>
<p><table cellpadding=10>

<tr><td>No parameter definitions for midplane_power.py.run_midplane_power

</table>

<h2>/home/phil/vita/core/vita/controller/midplane_to_divertor_power.py</h2>
<p><table cellpadding=10>

<tr><td>No parameter definitions for midplane_to_divertor_power.py.run_midplane_to_div_power

</table>

<h2>/home/phil/vita/core/vita/controller/parser.py</h2>
<p><table cellpadding=10>

<tr><td>No parameter definitions for parser.py.is_valid_file

</table>

<h2>/home/phil/vita/core/vita/modules/cherab/interface_surface.py</h2>
<p><table cellpadding=10>

<tr><td>No parameter definitions for InterfaceSurface.__init__

<tr><td>No parameter definitions for InterfaceSurface._generate_sample_point

<tr><td>No parameter definitions for InterfaceSurface._write_mesh_power_vtk

<tr><td>No parameter definitions for InterfaceSurface._write_mesh_points_vtk

<tr><td>Differing numbers of parameters described in comment and code

<tr><td>Parameter self not described by :param

<tr><td>Parameter equilibrium not described by :param

<tr><td>Parameters EFITEquilibrium equilibrium defined by :param but not present in defn

<tr><td>No parameter definitions for InterfaceSurface.plot

<tr><td>No parameter definitions for InterfaceSurface.histogram_plot

<tr><td>No parameter definitions for InterfaceSurface.get_lcfs_radius

<tr><td>No parameter definitions for InterfaceSurface.psin

</table>

<h2>/home/phil/vita/core/vita/modules/cherab/machine_configuration.py</h2>
<p><table cellpadding=10>

<tr><td>No parameter definitions for machine_configuration.py.load_wall_configuration

</table>

<h2>/home/phil/vita/core/vita/modules/equilibrium/equdsk/equdsk.py</h2>
<p><table cellpadding=10>

<tr><td>No parameter definitions for equdsk.py.read_eqdsk

</table>

<h2>/home/phil/vita/core/vita/modules/equilibrium/fiesta/fiesta_interface.py</h2>
<p><table cellpadding=10>

<tr><td>No parameter definitions for Fiesta().__init__

<tr><td>No parameter definitions for Fiesta().b_field

<tr><td>No parameter definitions for Fiesta().read_fiesta_model

<tr><td>No parameter definitions for Fiesta().get_midplane_lcfs

<tr><td>No parameter definitions for Fiesta().to_cherab_equilibrium

</table>

<h2>/home/phil/vita/core/vita/modules/equilibrium/fiesta/fiesta_interface_simple.py</h2>
<p><table cellpadding=10>

<tr><td>No parameter definitions for Fiesta().__init__

<tr><td>No parameter definitions for Fiesta().b_field

<tr><td>No parameter definitions for Fiesta().read_fiesta_model

<tr><td>No parameter definitions for Fiesta().get_midplane_lcfs

</table>

<h2>/home/phil/vita/core/vita/modules/projection/projection2D/field_line/field_line.py</h2>
<p><table cellpadding=10>

<tr><td>No parameter definitions for FieldLine().__init__

<tr><td>No parameter definitions for FieldLine().dx_dl

<tr><td>No parameter definitions for FieldLine().event

</table>

<h2>/home/phil/vita/core/vita/modules/projection/projection2D/field_line/field_line_projection.py</h2>
<p><table cellpadding=10>

<tr><td>No parameter definitions for field_line_projection.py._flux_expansion

<tr><td>No parameter definitions for field_line_projection.py.project_field_lines

</table>

<h2>/home/phil/vita/core/vita/modules/projection/projection2D/field_line/map_field_lines.py</h2>
<p><table cellpadding=10>

<tr><td>No parameter definitions for map_field_lines.py.map_field_lines

</table>

<h2>/home/phil/vita/core/vita/modules/projection/projection2D/particle_path_projection.py</h2>
<p><table cellpadding=10>

<tr><td>No parameter definitions for ParticlePath().__init__

<tr><td>No parameter definitions for ParticlePath().exact_b_phi

<tr><td>No parameter definitions for ParticlePath().dlorentz_dt

<tr><td>No parameter definitions for ParticlePath().dvec_dt

<tr><td>No parameter definitions for ParticlePath().event

<tr><td>No parameter definitions for ParticlePath().event2

</table>

<h2>/home/phil/vita/core/vita/modules/projection/projection2D/project_heat_flux.py</h2>
<p><table cellpadding=10>

<tr><td>No parameter definitions for project_heat_flux.py.project_heat_flux

</table>

<h2>/home/phil/vita/core/vita/modules/projection/projection2D/psi_map_projection.py</h2>
<p><table cellpadding=10>

<tr><td>No parameter definitions for psi_map_projection.py._flux_expansion

<tr><td>No parameter definitions for psi_map_projection.py._calculate_angles

<tr><td>No parameter definitions for psi_map_projection.py._angle

<tr><td>No parameter definitions for psi_map_projection.py.map_psi_omp_to_divertor

</table>

<h2>/home/phil/vita/core/vita/modules/sol_heat_flux/eich/eich.py</h2>
<p><table cellpadding=10>

<tr><td>No parameter definitions for Eich(HeatLoad).__init__

<tr><td>No parameter definitions for Eich(HeatLoad).calculate_heat_flux_density

</table>

<h2>/home/phil/vita/core/vita/modules/sol_heat_flux/hesel/hesel_data.py</h2>
<p><table cellpadding=10>

<tr><td>No parameter definitions for HESELdata().__init__

<tr><td>No parameter definitions for HESELdata().__del__

<tr><td>No parameter definitions for HESELdata()._load_file

<tr><td>No parameter definitions for HESELdata()._close_file

<tr><td>No parameter definitions for HESELdata()._evaluate_electron_advection_field

<tr><td>No parameter definitions for HESELdata()._evaluate_ion_advection_field

<tr><td>No parameter definitions for HESELdata()._evaluate_electron_conduction_field

<tr><td>No parameter definitions for HESELdata()._evaluate_ion_conduction_field

<tr><td>No parameter definitions for HESELdata()._evaluate_parallel_heat_flux_electron_advection

<tr><td>No parameter definitions for HESELdata()._evaluate_parallel_heat_flux_ion_advection

<tr><td>No parameter definitions for HESELdata()._evaluate_parallel_heat_flux_electron_conduction

<tr><td>No parameter definitions for HESELdata()._evaluate_parallel_heat_flux_ion_conduction

<tr><td>No parameter definitions for HESELdata().evaluate_parallel_heat_fluxes

<tr><td>No parameter definitions for HESELdata().get_lcfs_values

<tr><td>No parameter definitions for HESELdata().calculate_lambda_q

<tr><td>No parameter definitions for HESELdata().get_probe_positions

<tr><td>No parameter definitions for HESELdata()._load_probe_data

<tr><td>No parameter definitions for HESELdata().get_profiles_from_probes

<tr><td>No parameter definitions for HESELdata().load_2d_animation_fields

<tr><td>No parameter definitions for HESELdata().animate_2d_field

<tr><td>No parameter definitions for HESELdata().animate_1d_field

</table>

<h2>/home/phil/vita/core/vita/modules/sol_heat_flux/hesel/hesel_parameters.py</h2>
<p><table cellpadding=10>

<tr><td>No parameter definitions for HESELparams().__init__

<tr><td>No parameter definitions for HESELparams().read_hesel_output

<tr><td>No parameter definitions for HESELparams().calculate_sound_speed

<tr><td>No parameter definitions for HESELparams().calculate_ion_gyrofrequency

<tr><td>No parameter definitions for HESELparams().calculate_ion_gyroradius

<tr><td>No parameter definitions for HESELparams().calculate_debye_length

<tr><td>No parameter definitions for HESELparams().calculate_coulomb_logarithm

<tr><td>No parameter definitions for HESELparams().get_x_axis_rhos

<tr><td>No parameter definitions for HESELparams().get_y_axis_rhos

<tr><td>No parameter definitions for HESELparams().get_lcfs_index

<tr><td>No parameter definitions for HESELparams().get_wall_index

<tr><td>No parameter definitions for HESELparams().get_x_axis_probes

</table>

<h2>/home/phil/vita/core/vita/modules/sol_heat_flux/mid_plane_heat_flux.py</h2>
<p><table cellpadding=10>

<tr><td>No parameter definitions for HeatLoad().__init__

<tr><td>No parameter definitions for HeatLoad().__call__

<tr><td>No parameter definitions for HeatLoad().calculate_heat_flux_density

<tr><td>No parameter definitions for HeatLoad().set_coordinates

<tr><td>No parameter definitions for HeatLoad().set_edge_radius

<tr><td>No parameter definitions for HeatLoad().get_local_coordinates

<tr><td>No parameter definitions for HeatLoad().get_global_coordinates

<tr><td>No parameter definitions for HeatLoad().calculate_heat_power

<tr><td>No parameter definitions for HeatLoad().plot_heat_power_density

</table>

<h2>/home/phil/vita/core/vita/modules/utils/Vector2.py</h2>
<p><table cellpadding=10>

<tr><td>Number of parameters specified: 1 does not equal number of parameters documented: 1 on line: 23

<tr><td>Number of parameters specified: 0 does not equal number of parameters documented: 0 on line: 84

</table>

<h2>/home/phil/vita/core/vita/modules/utils/calculate_angle.py</h2>
<p><table cellpadding=10>

<tr><td>No parameter definitions for calculate_angle.py.calculate_angle

<tr><td>No parameter definitions for calculate_angle.py._angle

</table>

<h2>/home/phil/vita/core/vita/modules/utils/getOption.py</h2>
<p><table cellpadding=10>

<tr><td>No parameter definitions for getOption.py.getOption

</table>

<h2>/home/phil/vita/core/vita/modules/utils/intersections.py</h2>
<p><table cellpadding=10>

<tr><td>No parameter definitions for intersections.py._get_rectangle_intersections

<tr><td>No parameter definitions for intersections.py.intersection

</table>

<h2>/home/phil/vita/core/vita/modules/utils/load_pickle.py</h2>
<p><table cellpadding=10>

<tr><td>No parameter definitions for load_pickle.py.load_pickle

</table>

<h2>/home/phil/vita/core/vita/modules/utils/physics_constants.py</h2>
<p><table cellpadding=10>

<tr><td>No parameter definitions for physics_constants.py.get_physics_constants

</table>

<h2>/home/phil/vita/core/vita/utility/resource_manager.py</h2>
<p><table cellpadding=10>

<tr><td>No parameter definitions for resource_manager.py._test_allowed_characters

<tr><td>Parameter symlink=False not described by :param

<tr><td>Parameters symlink defined by :param but not present in defn

<tr><td>No parameter definitions for resource_manager.py.add_resource

<tr><td>Parameter symlink=False not described by :param

<tr><td>Parameters symlink defined by :param but not present in defn

<tr><td>No parameter definitions for resource_manager.py.update_resource

<tr><td>No parameter definitions for resource_manager.py.list_resources

<tr><td>No parameter definitions for resource_manager.py.get_resource

<tr><td>Parameter prompt_user=False not described by :param

<tr><td>Parameters prompt_user defined by :param but not present in defn

<tr><td>No parameter definitions for resource_manager.py.remove_resource

</table>

<h2>/home/phil/vita/core/vita/utility/scripts/add_resource.py</h2>
<p><table cellpadding=10>

<tr><td>No parameter definitions for add_resource.py.main

</table>

<h2>/home/phil/vita/core/vita/utility/scripts/remove_resource.py</h2>
<p><table cellpadding=10>

<tr><td>No parameter definitions for remove_resource.py.main

</table>

<h2>/home/phil/vita/core/vita/utility/scripts/update_resource.py</h2>
<p><table cellpadding=10>

<tr><td>No parameter definitions for update_resource.py.main

</table>

<h2>/home/phil/vita/core/vita/vitarunner.py</h2>
<p><table cellpadding=10>

<tr><td>No parameter definitions for vitarunner.py.main

</table>
