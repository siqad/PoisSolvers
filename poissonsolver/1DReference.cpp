//Reference poisson functions
void poisson1DJacobi(void) {
    //Parameters
    unsigned int N = 1000;           //Number of lattice points
    double sampleLength = 1;   //Physical length of sample in meters
    std::vector<double> V(N,0);     //Electric potential solution to Poisson's equation (Volts)
    std::vector<double> Vold(N,0);     //Vector to copy old values
    std::vector<double> rho(N,1.6e-19);   //Charge density, first and last are boundary conditions
    double maxError = 1e-5;       //set at 1% error allowed
    double currError;
    unsigned int cycleCount = 0;
    //Setup boundary condition
    std::cout << "Setup for Jacobi" << std::endl;
    rho[0] = 0;
    rho[N-1] = 0;
    //Prepare rho with EPS and spacing multiplication
    std::transform( rho.begin()+1, rho.end()-1,
                    rho.begin()+1, std::bind2nd(std::divides<double>(), EPS0 ) );
    std::transform( rho.begin()+1, rho.end()-1,
                    rho.begin()+1, std::bind2nd(std::multiplies<double>(),(sampleLength/(N-1)*(sampleLength/(N-1)))));
    //Apply coundary condition
    std::cout << "Applying boundary condition" << std::endl;
    V[0] = rho[0];
    V[N-1] = rho[N-1];
    //Computation loop
    std::cout << "Iterating..." << std::endl;
    do{
        currError = 0;  //reset error for every run
        Vold = V; //need to copy over all old values for Jacobi
        for( int i = 1; i < N-1; i++){ //for all lattice points except endpoints
            V[i] = (Vold[i-1] + Vold[i+1] + rho[i])/2.0; //calculate new potential
            if ( fabs((V[i] - Vold[i])/V[i]) > currError ){ //capture worst case error
                currError = fabs((V[i] - Vold[i])/V[i]);
            }
        }
        cycleCount++;
    }while( currError > maxError );
    std::cout << "Finished in " << cycleCount << " iterations within " <<
              maxError*100 << "% error." << std::endl;
}

void poisson1DGaussSeidel(void) {
    //Parameters
    unsigned int N = 1000;           //Number of lattice points
    double sampleLength = 1;       //Physical length of sample in meters
    std::vector<double> V(N,0);    //Electric potential solution to Poisson's equation (Volts)
    std::vector<double> rho(N,1.6e-19);   //Charge density, first and last are boundary conditions
    double Vold;                   //Vold only requires one element to calculate error
    double maxError = 1e-5;        //set at 1% error allowed
    double currError;
    unsigned int cycleCount = 0;
    //Setup
    std::cout << "Setup for Gauss-Seidel" << std::endl;
    rho[0] = 0;
    rho[N-1] = 0;
    std::transform( rho.begin()+1, rho.end()-1,
                    rho.begin()+1, std::bind2nd(std::divides<double>(), EPS0 ) );
    std::transform( rho.begin()+1, rho.end()-1, //multiply by spacing squared
                    rho.begin()+1, std::bind2nd(std::multiplies<double>(),(sampleLength/(N-1)*(sampleLength/(N-1)))));
    V[0] = rho[0];
    V[N-1] = rho[N-1];
    //Computation loop
    std::cout << "Iterating..." << std::endl;
    do{
        currError = 0;  //reset error for every run
        for( int i = 1; i < N-1; i=i+1 ){ //for all lattice points except endpoints,
            Vold = V[i]; //Save for error comparison
            V[i] = (V[i-1] + V[i+1] + rho[i])/2.0; //calculate new potential
            if ( fabs((V[i] - Vold)/V[i]) > currError ){ //capture worst case error
                currError = fabs((V[i] - Vold)/V[i]);
            }
        }
        cycleCount++;
    }while( currError > maxError );
    std::cout << "Finished in " << cycleCount << " iterations within " <<
              maxError*100 << "% error." << std::endl;
}

void poisson1DSOR(void) {
    //Parameters
    unsigned int N = 1000;           //Number of lattice points
    double sampleLength = 1;       //Physical length of sample in meters
    std::vector<double> V(N,0);    //Electric potential solution to Poisson's equation (Volts)
    std::vector<double> rho(N,1.6e-19);   //Charge density, first and last are boundary conditions
    double Vold;                   //Vold only requires one element to calculate error
    double maxError = 1e-5;        //set at 1% error allowed
    double overrelax = 1.85;        //overrelaxation parameter for SOR method
    double currError;
    unsigned int cycleCount = 0;
    //Setup
    std::cout << "Setup for SOR" << std::endl;
    rho[0] = 0;
    rho[N-1] = 0;
    std::transform( rho.begin()+1, rho.end()-1,
                    rho.begin()+1, std::bind2nd(std::divides<double>(), EPS0 ) );
    std::transform( rho.begin()+1, rho.end()-1, //multiply by spacing squared
                    rho.begin()+1, std::bind2nd(std::multiplies<double>(),(sampleLength/(N-1)*(sampleLength/(N-1)))));
    V[0] = rho[0];
    V[N-1] = rho[N-1];
    //Computation loop
    std::cout << "Iterating..." << std::endl;
    do{
        currError = 0;  //reset error for every run
        for( int i = 1; i < N-1; i=i+1 ){ //for all lattice points except endpoints,
            Vold = V[i]; //Save for error comparison
            V[i] = (1-overrelax)*V[i] + overrelax*(V[i-1] + V[i+1] + rho[i])/2.0; //calculate new potential
            if ( fabs((V[i] - Vold)/V[i]) > currError ){ //capture worst case error
                currError = fabs((V[i] - Vold)/V[i]);
            }
        }
        cycleCount++;
    }while( currError > maxError );
    std::cout << "Finished in " << cycleCount << " iterations within " <<
              maxError*100 << "% error." << std::endl;
}
