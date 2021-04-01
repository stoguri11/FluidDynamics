
import java.awt.* ;
import javax.swing.* ;

public class LBM {

    //final static int NITER = 30000 ;
    final static int NITER = 5000 ;

    final static int NX = 1080, NY = 720 ;  // Lattice dimensions
    final static int Q = 9 ;  // num states
     
    final static double uLB = 0.06 ;  // Inlet velocity in lattice units

    final static double Re = 330.0 ;  // Reynolds number

    // useful constants for weights

    final static double W0 = 4.0 / 9 ;
    final static double W1 = 1.0 / 9 ;
    final static double W2 = 1.0 / 36 ;

    final static int CELL_SIZE = 2 ; 
    final static int OUTPUT_FREQ = 100 ;

    static Display display = new Display() ;

    static int [] [] c = new int [Q] [2] ;  // Lattice velocities
    static double [] w = new double [Q] ;  // Lattice weights

    static double [] [] [] u = new double [NX] [NY] [2] ;

    static boolean [] [] obstacle = new boolean [NX] [NY] ;

    public static void main(String args []) throws Exception {

        double cx = NX/4.0, cy = NY/2.0, r = 20 ;
                // Coordinates and size of obstacle.

        double nulb = uLB * r / Re ;
        double omega = 1.0 / (3 * nulb + 0.5) ;  // Relaxation parameter

        long startTime = System.currentTimeMillis();

        long macroTime = 0 ;
        long collisionTime = 0 ;
        long streamingTime = 0 ;

        // Define table c of velocity states
        int [] vStates = new int [] {0, -1, +1} ;
        int pos = 0 ;
        for(int i : vStates)
            for(int j : vStates) {
                int [] cEl = c [pos++] ;
                cEl [0] = i ;
                cEl [1] = j ;
            }

        // Define table w of weights for equilibrium distribution
        w [0] = W0 ;
        for(int i = 1 ; i < Q ; i++) {
            int [] cEl = c [i] ;
            if(cEl [0] == 0 || cEl [1] == 0) {
                w [i] = W1 ;
            }
            else {
                w [i] = W2 ;
            }
        }

        int [] noslip = new int [Q] ;  // index in c of negative velocity state
        for(int i = 0 ; i < Q ; i++) {
            int [] cEl = c [i] ;
            for(int j = 0 ; j < Q ; j++) {
                int [] cElj = c [j] ;
                if(cElj [0] == -cEl [0] && cElj [1] == -cEl [1]) {
                    noslip [i] = j ;
                }
            }
        }

        int [] i1 = new int [3], i2 = new int [3], i3 = new int [3] ;
        int i1pos = 0, i2pos = 0, i3pos = 0 ;
        for(int i = 0 ; i < Q ; i++) {
            int cElX = c [i] [0] ;
            if (cElX < 0) {
                i1 [i1pos++] = i ; 
            }
            else if (cElX == 0) {
                i2 [i2pos++] = i ; 
            }
            else {
                i3 [i3pos++] = i ; 
            }
        }

        // Cylindrical obstacle
        double r2 = r * r ;
        for(int i = 0 ; i < NX ; i++) {
            for(int j = 0 ; j < NY ; j++) {
                double dx = i - cx ;
                double dy = j - cy ;
                obstacle [i] [j] = (dx * dx + dy * dy) < r2 ;
            }
        }

        // Inlet velocity with perturbation
        double [] [] vel = new double [NY] [2] ;
        for(int j = 0 ; j < NY ; j++) {
            vel [j] [0] = uLB *
                    (1 + 1E-4 * Math.sin(2 * Math.PI * j / (NY - 1))) ;
        }

        double [] [] [] fin = new double [NX] [NY] [Q] ;
        for(int i = 0 ; i < NX ; i++) {
            for(int j = 0 ; j < NY ; j++) {
                if(i == 0) {
                    equilibrium(fin [i] [j], 1.0, vel [j] [0], vel [j] [1]) ;
                }
                else {
                    equilibrium(fin [i] [j], 1.0, 0.0, 0.0) ;
                } 
            }
        }

        double [] [] [] fout = new double [NX] [NY] [Q] ;
        double [] [] rho = new double [NX] [NY] ;

        for(int time = 0 ; time < NITER ; time++) {

            long time1 = System.currentTimeMillis();

            // Calculate macroscopic density and velocity
            for(int i = 0 ; i < NX ; i++) {
                for(int j = 0 ; j < NY ; j++) {
                    double [] fin_ij = fin [i] [j] ;
                    double [] u_ij = u [i] [j] ;
                    if(i > 0) {
                        double sum = fin_ij [0] + fin_ij [1] + fin_ij [2] +
                                     fin_ij [3] + fin_ij [4] + fin_ij [5] +
                                     fin_ij [6] + fin_ij [7] + fin_ij [8] ;

                        double sum0 = - fin_ij [3] - fin_ij [4] - fin_ij [5]
                                      + fin_ij [6] + fin_ij [7] + fin_ij [8] ;

                        double sum1 = - fin_ij [1] + fin_ij [2] - fin_ij [4]
                                      + fin_ij [5] - fin_ij [7] + fin_ij [8] ;


                        rho [i] [j] = sum ;
                        if(sum > 0) {
                            u_ij [0] = sum0 / sum ;
                            u_ij [1] = sum1 / sum ;
                        }
                    }
                    else {
                        // BC - left wall: compute density from known
                        // populations.

                        u_ij [0] = vel [j] [0] ;
                        u_ij [1] = vel [j] [1] ;
                        float sum2 = 0 ;
                        for(int d : i2) {
                           sum2 += fin_ij [d] ;
                        }
                        float sum1 = 0 ;
                        for(int d : i1) {
                           sum1 += fin_ij [d] ;
                        }
                        rho [0] [j] = 1/(1 - u_ij [0]) * (sum2 + 2 * sum1) ;
                    }
                }
            }

            long time2 = System.currentTimeMillis();

            macroTime += (time2 - time1) ;

            // Collision step.
            for(int i = 0 ; i < NX ; i++) {
                for(int j = 0 ; j < NY ; j++) {
                    double [] fin_ij = fin [i] [j] ;
                    double [] fout_ij = fout [i] [j] ;
                    if(obstacle [i] [j]) {
                        // BC - no slip at obstacle
                        for(int d = 0; d < Q ; d++) {
                            fout_ij [d] = fin_ij [noslip [d]] ;
                        }
                    }
                    else {
                        double [] feq = new double [Q] ;

                        equilibrium(feq,
                                    rho [i] [j], u [i] [j] [0], u [i] [j] [1]) ;

                        // BC - Left wall: Equilibrium scheme
                        if(i == 0) {
                            for(int p = 0 ; p < 3 ; p++) {
                                fin_ij [i3 [p]] = feq [i3 [p]] ;
                            }
                        }
                        fout_ij [0] = fin_ij [0] -
                                      omega * (fin_ij [0] - feq [0]) ;
                        fout_ij [1] = fin_ij [1] -
                                      omega * (fin_ij [1] - feq [1]) ;
                        fout_ij [2] = fin_ij [2] -
                                      omega * (fin_ij [2] - feq [2]) ;
                        fout_ij [3] = fin_ij [3] -
                                      omega * (fin_ij [3] - feq [3]) ;
                        fout_ij [4] = fin_ij [4] -
                                      omega * (fin_ij [4] - feq [4]) ;
                        fout_ij [5] = fin_ij [5] -
                                      omega * (fin_ij [5] - feq [5]) ;
                        fout_ij [6] = fin_ij [6] -
                                      omega * (fin_ij [6] - feq [6]) ;
                        fout_ij [7] = fin_ij [7] -
                                      omega * (fin_ij [7] - feq [7]) ;
                        fout_ij [8] = fin_ij [8] -
                                      omega * (fin_ij [8] - feq [8]) ;
                    }
                }
            }

            long time3 = System.currentTimeMillis();

            collisionTime += (time3 - time2) ;

            // Streaming step.
            for(int i = 0 ; i < NX ; i++) {

                int iP1 = (i + 1) % NX ;
                int iM1 = (i - 1 + NX) % NX ;

                double [] [] fin_i = fin [i] ;
                double [] [] fin_iM1 = fin [iM1] ;
                double [] [] fin_iP1 = fin [iP1] ;

                for(int j = 0 ; j < NY ; j++) {
                    double [] fout_ij = fout [i] [j] ;
                    int jP1 = (j + 1) % NY ;
                    int jM1 = (j - 1 + NY) % NY ;

                    fin_i [j] [0] = fout_ij [0] ;
                    fin_i [jM1] [1] = fout_ij [1] ;
                    fin_i [jP1] [2] = fout_ij [2] ;
                    fin_iM1 [j] [3] = fout_ij [3] ;
                    fin_iM1 [jM1] [4] = fout_ij [4] ;
                    fin_iM1 [jP1] [5] = fout_ij [5] ;
                    fin_iP1 [j] [6] = fout_ij [6] ;
                    fin_iP1 [jM1] [7] = fout_ij [7] ;
                    fin_iP1 [jP1] [8] = fout_ij [8] ;
                }
            }

            // BC - Right wall: outflow condition
            for(int d : i1) {
                for(int j = 0 ; j < NY ; j++) {
                    fin [NX - 1] [j] [d] = fin [NX - 2] [j] [d] ;
                }
            }

            long time4 = System.currentTimeMillis();

            streamingTime += (time4 - time3) ;

            if(time % OUTPUT_FREQ == 0) {
                System.out.println("time = " + time + "/" + NITER) ;
                display.repaint() ;
            }
        }


        long endTime = System.currentTimeMillis();

        System.out.println("Calculation completed in " +
                           (endTime - startTime) + " milliseconds");

        System.out.println("Time to calculate macroscopic quantities: " +
                           macroTime + " milliseconds");
        System.out.println("Time for collision steps: " +
                           collisionTime + " milliseconds");
        System.out.println("Time for streaming steps: " +
                           streamingTime + " milliseconds");

        display.repaint() ;
    }

    
    static void equilibrium(double [] feq, double rho, double u0, double u1) {

        double usqr = u0 * u0 + u1 * u1 ;

        // UNROLLED version of above loop over d
        double u0Pu1 = u0 + u1 ;
        double u0Mu1 = u0 - u1 ;
 
        feq [0] = rho * W0 * (1.0 - 1.5 * usqr) ;
        feq [1] = rho * W1 * (1.0 - 3.0 * u1 + 4.5 * u1 * u1 - 1.5 * usqr) ;
        feq [2] = rho * W1 * (1.0 + 3.0 * u1 + 4.5 * u1 * u1 - 1.5 * usqr) ;
        feq [3] = rho * W1 * (1.0 - 3.0 * u0 + 4.5 * u0 * u0 - 1.5 * usqr) ;
        feq [4] = rho * W2 * (1.0 - 3.0 * u0Pu1 +
                              4.5 * u0Pu1 * u0Pu1 - 1.5 * usqr) ;
        feq [5] = rho * W2 * (1.0 - 3.0 * u0Mu1 +
                              4.5 * u0Mu1 * u0Mu1 - 1.5 * usqr) ;
        feq [6] = rho * W1 * (1.0 + 3.0 * u0 + 4.5 * u0 * u0 - 1.5 * usqr) ;
        feq [7] = rho * W2 * (1.0 + 3.0 * u0Mu1 +
                              4.5 * u0Mu1 * u0Mu1 - 1.5 * usqr) ;
        feq [8] = rho * W2 * (1.0 + 3.0 * u0Pu1 +
                              4.5 * u0Pu1 * u0Pu1 - 1.5 * usqr) ;
    }

    static class Display extends JPanel {

        Display() {

            setPreferredSize(new Dimension(CELL_SIZE * NX, CELL_SIZE * NY)) ;

            JFrame frame = new JFrame("LBM");
            frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
            frame.setContentPane(this);
            frame.pack();
            frame.setVisible(true);
        }

        public void paintComponent(Graphics g) {
            double [] [] usqr = new double [NX] [NY] ;
            double uMax = Double.MIN_VALUE ;
            double uMin = Double.MAX_VALUE ;
            for(int i = 0 ; i < NX ; i++) {
                for(int j = 0 ; j < NY ; j++) {
                   double u0 = u [i] [j] [0] ;
                   double u1 = u [i] [j] [1] ;
                   double u2 = Math.sqrt(u0 * u0 + u1 * u1) ;
                   if(u2 < uMin) {
                       uMin = u2 ;
                   } 
                   if(u2 > uMax) {
                       uMax = u2 ;
                   } 
                   usqr [i] [j] = u2 ;
               
                }
            }
            double diff = uMax - uMin ;
            double norm = ((diff == 0.0) ? 0.0 : 1/diff) ;
            for(int i = 0 ; i < NX ; i++) {
                for(int j = 0 ; j < NY ; j++) {
                    if(obstacle [i] [j]) {
                        g.setColor(Color.BLACK) ;
                    }
                    else {
                        float f = (float) (norm * (usqr [i] [j] - uMin)) ;
                        Color c = new Color(f, 0.0F, 1.0F - f) ;
                        g.setColor(c) ;
                    }
                    g.fillRect(CELL_SIZE * i, CELL_SIZE * j,
                               CELL_SIZE, CELL_SIZE) ;
                }
            } 
            g.setColor(Color.WHITE) ;
            for(int i = 0 ; i < NX ; i+=8) {
                for(int j = 0 ; j < NY ; j+=8) {
                    int originX = CELL_SIZE * i ;
                    int originY = CELL_SIZE * j ;
                    g.drawOval(originX - 1, originY - 1, 3, 3) ;
                    g.drawLine(originX, originY,
                               originX + (int) (200 * u [i] [j] [0]),
                               originY + (int) (200 * u [i] [j] [1])) ;
                }
            } 
        }
    }
}
