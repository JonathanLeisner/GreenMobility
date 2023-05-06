""" 
Defines the class that runs the Green Mobility (GM) project code. 
The straightforward way of running the tests is to open a console and write
'python -m unittest -v tests'
or
'python -m unittest -v tests.TestGradients_partial' 
to run a specific set of tests
or
'python -m unittest -v tests.TestGradients_partial.test_estimation_from_shuffled_values_all_parameters' 
to run a specific individual test.
"""
#%%
import unittest
import numpy as np
from scipy import optimize
import util
from estimate import GM_estimate
#%%

class TestGradients_partial(unittest.TestCase):
    """
    Class for testing gradients and correct estimation in the partial version of the model, 
    i.e. where skill prices (r) are constant throughout.
    """
    @classmethod
    def setUpClass(cls):
        """
        Prepare a class instance of a partial version of the model, as a pickle"""
        gm = GM_estimate(name="TestGradients_partial", endo_D=True)
        gm.partial_model = True
        gm.setup()
        gm.allocate_worker()
        gm.precompute()
        gm.solve_worker()
        gm.setup_simulated_data(n_individuals=100_000)
        gm.to_pickle()

    def setUp(self):
        """ 
        Loads a model instance from a pickle and sets the standard error tolerance for evaluating whether
        estimation comes sufficiently close to the true values.
        """
        self.gm = GM_estimate(from_pickle="TestGradients_partial", endo_D=False)
        self.SE_tol = 3

    def assertGradients(self):
        """ 
        If level difference is small enough, accept. If not, check if relative difference is small enough.
        """
        df = self.gm.compare_gradients_loglik()
        msg = f"\nGradient comparison: \n {df}"
        maxdiff = np.max(np.abs(df["Difference"]))
        if maxdiff <= 1e-6:
            self.assertLess(maxdiff, 1e-6, msg)
        else:
            self.assertLess(np.max(np.abs(df["Relative difference"])), 0.05, msg)

    def test_analytic_vs_numeric_gradient_individual_parameters(self):
        """
        Tests whether the analytic and numeric gradients for individual parameters are sufficiently close.
        Note that 'individual parameters' refers to individual names, i.e. beta0. So this tests the 
        S - 1 derivatives in beta0, corresponding to each of the elements in beta0.
        """
        gm = self.gm
        for p in gm.par.all_parameters:
            gm.reset_par()
            gm.setup_estimation(parameters=[p], agrad_loglik=False)
            with self.subTest(p=p):
                self.assertGradients()

    def test_analytic_vs_numeric_gradient_all_parameters(self):
        """
        Tests whether all analytic and numeric gradients are sufficiently close when calculated jointly for all parameters.
        """
        gm = self.gm
        gm.reset_par()
        gm.setup_estimation(agrad_loglik=False)
        self.assertGradients()

    def test_gradient_increases_away_from_true_values(self):
        """
        Tests that the analytic gradient increases when the point that the derivative is taken in is moved away from the true optimum.
        """
        gm = self.gm
        for p in gm.par.all_parameters:
            for idx in np.arange(getattr(gm.par, p).size):
                if p in gm.par.restrictions:
                    if idx in gm.par.restrictions[p]:
                        continue
                gm.reset_par()
                gm.setup_estimation(parameters=[p], indexes=np.array([idx]), agrad_loglik=False)
                theta0 = gm.est.theta0.copy()
                grad_in_optimum = gm.score_from_theta(theta0)
                for pertubation in [0.4, 1.6]:
                    gm.reset_par()
                    with self.subTest(parameter=p, pertubation=pertubation):
                        diff = np.abs(gm.score_from_theta(theta0 * pertubation)) - np.abs(grad_in_optimum)
                        self.assertGreater(diff, 0)

    def assertEstimation(self):
        """ 
        Checks whether an estimation results in parameter estimates that are sufficiently close to the true ones. 
        """
        df = self.gm.compare_init_to_estimated()
        msg = f"\nPost estimation summary: \n{df}\n"
        self.assertLess(np.max(np.abs(df["t-stat (p=default)"])), self.SE_tol, msg)

    def test_estimation_from_true_values_individual_parameters(self):
        """
        Tests whether estimation results in parameter estimates that are sufficiently close to the true ones.
        The starting values are the true ones themselves. Tests individual parameters separately e.g. beta0.
        """
        gm = self.gm
        for p in gm.par.all_parameters:
            gm.reset_par()
            gm.setup_estimation(parameters=[p], agrad_loglik=True)
            gm.estimate()
            self.assertEstimation()

    def test_estimation_from_shuffled_values_individual_parameters(self):
        """
        Tests whether estimation results in parameter estimates that are sufficiently close to the true ones.
        The starting values are moved away from the true ones. Tests individual parameters separately e.g. beta0.
        """
        gm = self.gm
        for p in gm.par.all_parameters:
            gm.reset_par()
            gm.setup_estimation(parameters=[p], agrad_loglik=True)
            gm.est.theta0 = gm.est.theta0 * (util.shuffle_vector(len(gm.est.theta0)) + 1)
            gm.estimate()
            self.assertEstimation()

    def test_estimation_from_true_values_all_parameters(self):
        """
        Tests whether estimation results in parameter estimates that are sufficiently close to the true ones.
        The starting values are the true ones. Tests all parameters jointly.
        """
        gm = self.gm
        gm.reset_par()
        gm.setup_estimation(agrad_loglik=True)
        gm.estimate()
        self.assertEstimation()
        gm.t_compare_init_to_estimated()
        gm.tables_to_pdf(gm.name + "estimates_from_true_values")

    def test_estimation_from_shuffled_values_all_parameters(self):
        """
        Tests whether estimation results in parameter estimates that are sufficiently close to the true ones.
        The starting values are moved away from the true ones. Tests all parameters jointly.
        """
        gm = self.gm
        gm.reset_par()
        gm.setup_estimation(agrad_loglik=True)
        gm.est.theta0 = gm.est.theta0 * (util.shuffle_vector(len(gm.est.theta0)) + 1)
        gm.estimate()
        gm.t_compare_init_to_estimated()
        gm.tables_to_pdf(gm.name + "estimates_from_shuffled_values")
        self.assertEstimation()
        

class TestGradients_r(unittest.TestCase):
    """
    Class for testing gradients with respect to skill prices (r). These are used for reaching labor market
    equilibrium and in calculation of derivatives wrt. the likelihood function when the model is not in partial mode.
    """
    @classmethod
    def setUpClass(cls):
        """
        Prepares a pickle of a model instance that is not in equilibrium. 
        """
        gm = GM_estimate(name="TestGradients_r", endo_D=True)
        gm.partial_model = True
        gm.setup()
        gm.allocate_worker()
        gm.precompute()
        gm.solve_worker()
        gm.setup_simulated_data(n_individuals=100_000)
        gm.partial_model = False
        gm.to_pickle()

    def setUp(self):
        """
        Loads the model instance as preparation for running a test on it.
        """
        self.gm = GM_estimate(from_pickle="TestGradients_r", endo_D=False)

    def test_inner_equilibrium_gradients_outofequilibrium(self):
        """
        Tests whether the gradients of the excess labor demand functions with respect to prices r are equal, 
        whether calculated numerically or analytically. They are evaluated in the initial values, that is, out of equilibrium.
        """
        gm = self.gm
        r_1d = gm.par.r.reshape((gm.par.T * (gm.par.S - 1)), order="F")
        n = 0
        for s in np.arange(gm.par.S - 1):
            for t in np.arange(gm.par.T):
                f = lambda x: gm.ED_from_r(x)[n]
                g = lambda x: gm.dED_dr(x)[t, s, :, :].reshape((gm.par.T * (gm.par.S - 1)), order="F")
                num = optimize.approx_fprime(xk=r_1d, f=f, epsilon=1.4901161193847656e-08)
                ana = g(r_1d)
                with self.subTest(sector=s, year=t):
                    self.assertLess(np.max(np.abs(num - ana)), 1e-6)
                n += 1

    def test_inner_equilibrium_gradients_inequilibrium(self):
        """
        Tests whether the gradients of the excess labor demand functions with respect to prices r are equal, 
        whether calculated numerically or analytically. They are evaluated in their equilibrium values.
        """
        gm = self.gm
        gm.find_humcap_equilibrium()
        r_1d = gm.par.r.reshape((gm.par.T * (gm.par.S - 1)), order="F")
        n = 0
        for s in np.arange(gm.par.S - 1):
            for t in np.arange(gm.par.T):
                f = lambda x: gm.ED_from_r(x)[n]
                g = lambda x: gm.dED_dr(x)[t, s, :, :].reshape((gm.par.T * (gm.par.S - 1)), order="F")
                num = optimize.approx_fprime(xk=r_1d, f=f, epsilon=1.4901161193847656e-08)
                ana = g(r_1d)
                with self.subTest(sector=s, year=t):
                    self.assertLess(np.max(np.abs(num - ana)), 1e-6)
                n += 1

    def test_equilibrium_skill_prices_numeric_versus_analytic_solution(self):
        """
        Tests whether the equilibrium skill prices found without analytic gradients are the same as those found with
        numeric gradients.
        """
        gm = self.gm
        gm.agrad_quadED = True
        gm.agrad_ED0 = True
        gm.find_humcap_equilibrium()
        r_analytic = gm.par.r.copy()
        #Reset model:
        gm = GM_estimate(from_pickle="TestGradients_r", endo_D=False)
        gm.agrad_quadED = False
        gm.agrad_ED0 = False
        gm.find_humcap_equilibrium()
        r_numeric = gm.par.r.copy()
        self.assertLess(np.max(np.abs(r_analytic - r_numeric)), 1e-8)

class TestGradients_GE(unittest.TestCase):
    """
    Class for testing gradients and estimation results when skill prices are determined in equilibrium.
    """
    @classmethod
    def setUpClass(cls):
        """
        Prepares a pickle of a model instance that is in equilibrium. 
        """
        gm = GM_estimate(name="TestGradients_GE", endo_D=True)
        gm.setup()
        gm.allocate_worker()
        gm.precompute()
        gm.solve_worker()
        gm.setup_simulated_data(n_individuals=50_000)
        gm.simulate()
        #a lot of iterations to make sure we are very close to the stationary point
        gm.successive_approximations(95)
        gm.successive_approximations(95)
        gm.setup_simulated_data(n_individuals=100_000)
        gm.to_pickle()

    def setUp(self):
        """ 
        Loads a model instance where the model is in equilibrium and defines the standard error tolerance for evaluating
        estimation results.
        """ 
        self.gm = GM_estimate(from_pickle="TestGradients_GE", endo_D=False)
        self.SE_tol = 4

    def assertGradients(self):
        """ If level difference is small enough, accept. If not, check if relative difference is small enough."""
        df = self.gm.compare_gradients_loglik()
        msg = f"\nGradient comparison: \n {df}"
        maxdiff = np.max(np.abs(df["Difference"]))
        if maxdiff <= 1e-6:
            self.assertLess(maxdiff, 1e-6, msg)
        else:
            self.assertLess(np.max(np.abs(df["Relative difference"])), 0.05, msg)

    def test_analytic_vs_numeric_gradient_individual_parameters(self):
        """
        Tests whether the analytic and numeric gradients for individual parameters are sufficiently close.
        Note that 'individual parameters' refers to individual names, i.e. beta0. So this tests e.g. the 
        S - 1 derivatives in beta0, corresponding to each of the elements in beta0.
        """
        gm = self.gm
        for p in gm.par.all_parameters:
            gm.reset_par()
            gm.setup_estimation(parameters=[p], agrad_loglik=False)
            with self.subTest(p=p):
                self.assertGradients()

    def test_analytic_vs_numeric_gradient_all_parameters(self):
        """
        Tests whether the analytic and numeric gradients for all parameters jointly are sufficiently close.
        """
        gm = self.gm
        gm.reset_par()
        gm.setup_estimation(agrad_loglik=False)
        self.assertGradients()

    def test_gradient_increases_away_from_true_values(self):
        """ 
        Tests that the analytic gradient increases when the point that the derivative is taken in is moved away from the true optimum.
        """
        gm = self.gm
        for p in gm.par.all_parameters:
            for idx in np.arange(getattr(gm.par, p).size):
                if p in gm.par.restrictions:
                    if idx in gm.par.restrictions[p]:
                        continue
                gm.reset_par()
                gm.setup_estimation(parameters=[p], indexes=np.array([idx]), agrad_loglik=False)
                theta0 = gm.est.theta0.copy()
                grad_in_optimum = gm.score_from_theta(theta0)
                for pertubation in [0.4, 1.6]:
                    gm.reset_par()
                    with self.subTest(parameter=p, pertubation=pertubation):
                        diff = np.abs(gm.score_from_theta(theta0 * pertubation)) - np.abs(grad_in_optimum)
                        self.assertGreater(diff, 0)

    def assertEstimation(self):
        """
        Checks whether an estimation results in parameter estimates that are sufficiently close to the true ones. 
        """
        df = self.gm.compare_init_to_estimated()
        msg = f"\nPost estimation summary: \n{df}\n"
        self.assertLess(np.max(np.abs(df["t-stat (p=default)"])), self.SE_tol, msg)

    def test_estimation_from_true_values_individual_parameters(self):
        """
        Tests whether estimation results in parameter estimates that are sufficiently close to the true ones.
        The starting values are the true ones themselves. Tests individual parameters separately e.g. beta0.
        """
        gm = self.gm
        for p in gm.par.all_parameters:
            gm.reset_par()
            gm.setup_estimation(parameters=[p], agrad_loglik=True)
            gm.estimate()
            self.assertEstimation()

    def test_estimation_from_shuffled_values_individual_parameters(self):
        """
        Tests whether estimation results in parameter estimates that are sufficiently close to the true ones.
        The starting values are moved away from the true ones. Tests individual parameters separately e.g. beta0.
        """
        gm = self.gm
        for p in gm.par.all_parameters:
            gm.reset_par()
            gm.setup_estimation(parameters=[p], agrad_loglik=True)
            gm.est.theta0 = gm.est.theta0 * (util.shuffle_vector(len(gm.est.theta0)) + 1)
            gm.estimate()
            self.assertEstimation()

    def test_estimation_from_true_values_all_parameters(self):
        """
        Tests whether estimation results in parameter estimates that are sufficiently close to the true ones.
        The starting values are the true ones. Tests all parameters jointly.
        """
        gm = self.gm
        gm.reset_par()
        gm.setup_estimation(agrad_loglik=True)
        gm.estimate()
        gm.t_compare_init_to_estimated()
        gm.tables_to_pdf(gm.name + "estimates_from_true_values")
        self.assertEstimation()

    def test_estimation_from_shuffled_values_all_parameters(self):
        """
        Tests whether estimation results in parameter estimates that are sufficiently close to the true ones.
        The starting values are moved away from the true ones. Tests all parameters jointly.
        """
        gm = self.gm
        gm.reset_par()
        gm.setup_estimation(agrad_loglik=True)
        gm.est.theta0 = gm.est.theta0 * (util.shuffle_vector(len(gm.est.theta0)) + 1)
        gm.estimate()
        gm.t_compare_init_to_estimated()
        gm.tables_to_pdf(gm.name + "estimates_from_shuffled_values")
        self.assertEstimation()
