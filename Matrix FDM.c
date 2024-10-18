#include <stdio.h>
#include <stdlib.h>
#include <math.h>

//PART 1

//Global matrix containing deflections 
double **deflection;

//Global array containing locations of deflections
double *deflection_location;

//Local arrays 
double *deflection_left_edge;
double *deflection_right_edge;
double *deflection_upper_edge;
double *deflection_lower_edge; 

void initialize_edges(void);

double NDDM (double, double *); //This will be the boundary function f(a), with needed inputs: a and the specific edge
double divided_difference(int, int, double *);

void print_boundary_conditions(void);

//PART 2

double **final_deflection;

//For Gauss-Seidel Method
double **A;
double *B; 
double *x;

void transform_matrix(int);

double *gauss_seidel(double **A, double *B, int no_of_unknowns, double tolerance);

//PART 3

double *final_deflection_location;

//this function will retrieve the points of zero deflection given a deflection field
void get_points(int, double, double);

//these 2 arrays will contain the (x,y) coordinates of the points of zero deflection,
//while no_of_points gives the number of (x,y) points
int no_of_points;
double *x_zero;
double *y_zero;

//this function returns an array with three values, [0] = a, [1] = b, and [2] = r^2
//where the best fit line is y=ax+b, with goodness of fit r^2
double *least_squares(int N, double *x_values, double *y_values);

int main (void)
{
    //PART 1

    double mew = 5; //Average of last digit of our student numbers, rounded up to nearest int

    //Set empty matrix for deflection matrix
    //5 rows, 5 columns
    deflection = malloc(5*sizeof(double *));
    for (int i = 0; i < 5; i++)
    {
        deflection[i] = malloc(5*sizeof(double));
    }

    //Initialize values for deflection matrix    
    deflection[0][0] = 0; //corners
    deflection[0][4] = 0;
    deflection[4][0] = 0;
    deflection[4][4] = 0;

    deflection[0][1] = -0.092*mew; //left edge
    deflection[0][2] = -0.064*mew;
    deflection[0][3] = 0.011*mew;

    deflection[4][1] = -0.090*mew; //right edge
    deflection[4][2] = 0.055*mew;
    deflection[4][3] = 0.077*mew; 

    deflection[1][0] = -0.034*mew; //lower edge
    deflection[2][0] = 0.087*mew;
    deflection[3][0] = 0.046*mew;

    deflection[1][4] = -0.079*mew; //upper edge
    deflection[2][4] = -0.060*mew;
    deflection[3][4] = 0.090*mew;

    //Set empty matrix for deflection locations matrix
    deflection_location = malloc(5*sizeof(double));
    
    //Initialize values for deflection location matrix
    deflection_location[0] = 0;
    deflection_location[1] = 0.25*mew;
    deflection_location[2] = 0.5*mew;
    deflection_location[3] = 0.75*mew;
    deflection_location[4] = mew;

    initialize_edges();

    print_boundary_conditions();

    //PART 2

    //Ask for step size
    double step_size;
    printf("Enter step size for d_x and d_y (must divide perfectly into %.1lf): ", mew);
    scanf("%lf", &step_size);

    int nodes = (int) ((mew / step_size) + 1); 
    //if the deflection field has n^2 values, then there are n nodes

    puts("");
    printf("With a step size of %.3lf, this program will display %d (%d by %d) points.\n", step_size, (int) pow(nodes, 2), nodes, nodes);

    //Set empty matrix for final deflections matrix
    final_deflection = calloc(nodes, sizeof(double *));
    for (int i = 0; i < nodes; i++)
    {
        final_deflection[i] = calloc(nodes, sizeof(double));
    }

    //Using NDDM, we can now initialize values at corners and edges (using boundary conditions)
    final_deflection[0][0] = 0;
    final_deflection[0][nodes-1] = 0;
    final_deflection[nodes-1][0] = 0;
    final_deflection[nodes-1][nodes-1] = 0;

    for (int i = 1; i < nodes - 1; i++)
    {
        final_deflection[0][i] = NDDM((i/(double)(nodes-1)) * mew, deflection_left_edge); //left edge
        final_deflection[nodes-1][i] = NDDM((i/(double)(nodes-1)) * mew, deflection_right_edge); //right edge
        final_deflection[i][nodes-1] = NDDM((i/(double)(nodes-1)) * mew, deflection_upper_edge); //upper edge
        final_deflection[i][0] = NDDM((i/(double)(nodes-1)) * mew, deflection_lower_edge); //lower edge
    }

    printf("Using Gauss-Seidel Method, this program will solve for the inner %d (%d by %d) unknown points.\n", (int) pow(nodes - 2, 2), nodes - 2, nodes - 2);
    puts("");

    //Transform matrix - this will give us matrices A and B for gauss-seidel 
    transform_matrix(nodes);
        
    //Ask for tolerance for gauss-seidel method
    double tol;
    printf("Enter tolerance value: ", pow(nodes - 2, 2));
    scanf("%lf", &tol);
    puts("");

    //Solve unknowns through gauss-seidel method
    double *inner_unknowns = gauss_seidel(A, B, (int) pow(nodes - 2, 2), tol);

    //Transfer solved unknowns to final matrix
    int x_index = 0;

    for (int i = 1; i < (nodes-2) + 1; i++)
    {
        for (int j = 1; j < (nodes-2) + 1; j++)
        {
            final_deflection[j][i] = inner_unknowns[x_index];

            x_index++;
        }
    }

    //Print final matrix
    puts("");
    printf("Using Finite Difference Method, the steady state of the deflection field is as follows.\n");

    puts("");

    puts("y+");
    puts("^");
    puts("|");
    puts("|");
    puts("  ---> x+");
    puts("");

    for (int i = nodes-1; i > -1; i--)
    {
        for (int j = 0; j < nodes; j++)
        {
            printf("%.5lf \t", final_deflection[j][i]);
        }
        puts("");
        puts("");
    }

    //PART 3

    final_deflection_location = malloc(nodes * sizeof(double));
    for (int i = 0; i < nodes; i++)
    {
        final_deflection_location[i] = i * step_size;
    }

    get_points(nodes, step_size, mew);

    puts("The coordinates (x, y) of points of zero deflection are as follows:");
    printf("NOTE: These are limited to points such that %.3lf < x < %.3lf.\n", 0.3*mew, 0.7*mew);
    puts("");

    for (int i = 0; i < no_of_points; i++)
    {
        printf("(%.4lf, %.4lf)\n", x_zero[i], y_zero[i]);
    }

    double *best_fit_curve = least_squares(no_of_points, x_zero, y_zero);

    puts("");
    if (best_fit_curve[1] > 0)
    {
        printf("Using the least squares method, the equation for the best fit linear curve is y = %lfx + %lf.\nIt has an r^2 value of %lf.", best_fit_curve[0], best_fit_curve[1], best_fit_curve[2]);
    }
    else if (best_fit_curve[1] < 0)
    {
        printf("Using the least squares method, the equation for the best fit linear curve is y = %lfx - %lf.\nIt has an r^2 value of %lf.", best_fit_curve[0], fabs(best_fit_curve[1]), best_fit_curve[2]);
    }

    puts("");

    free(deflection);
    free(deflection_location);
    free(deflection_left_edge);
    free(deflection_right_edge);
    free(deflection_upper_edge);
    free(deflection_lower_edge);
    free(final_deflection);
    free(A);
    free(B);
    free(x);
    free(final_deflection_location);
    free(x_zero);
    free(y_zero);

    return 0;
}

void initialize_edges(void)
{
    deflection_left_edge = malloc(5*sizeof(double)); 
    for (int i = 0; i < 5; i++)
    {
        deflection_left_edge[i] = deflection[0][i];
    }

    deflection_right_edge = malloc(5*sizeof(double)); 
    for (int i = 0; i < 5; i++)
    {
        deflection_right_edge[i] = deflection[4][i];
    }

    deflection_upper_edge = malloc(5*sizeof(double)); 
    for (int i = 0; i < 5; i++)
    {
        deflection_upper_edge[i] = deflection[i][4];
    }

    deflection_lower_edge = malloc(5*sizeof(double)); 
    for (int i = 0; i < 5; i++)
    {
        deflection_lower_edge[i] = deflection[i][0];
    }

    return;
}

double NDDM (double x, double *I)
{
    double running_sum = 0;
    double running_termproduct = 1;

    for (int i = 0; i < 5; i++)
    {
        for (int j = 0; j < i; j++)
        {
            running_termproduct *= (x - deflection_location[j]);
        }

        running_sum += divided_difference(0, i, I) * running_termproduct;
         
        running_termproduct = 1;
    }

    return running_sum; 
}

double divided_difference(int i, int j, double* y_values)
{
    if (i == j) //base case
    {
        return y_values[i];
    }
    else
    {
        return (divided_difference(i+1, j, y_values) - divided_difference(i, j-1, y_values)) / (deflection_location[j] - deflection_location[i]);
    }
}

void print_boundary_conditions(void)
{
    printf("Using 4th Order Newton's Divided Difference Method, the boundary conditions for each edge may be approximated by: \n");
    puts("");

    printf("Left Edge: \n");

    printf("f(a) = ");

    for (int i = 0; i < 5; i++)
    {
        if (divided_difference(0, i, deflection_left_edge) > 0)
        {
            printf(" + %.3lf", divided_difference(0, i, deflection_left_edge));
        }
        else if (divided_difference(0, i, deflection_left_edge) < 0)
        {
            printf(" - %.3lf", fabs(divided_difference(0, i, deflection_left_edge)));
        }
        
        for (int j = 0; j < i; j++)
        {
            printf("(a-%.2lf)", deflection_location[j]);
        }
    }

    puts("");
    puts("");

    printf("Right Edge: \n");

    printf("f(a) = ");

    for (int i = 0; i < 5; i++)
    {
        if (divided_difference(0, i, deflection_right_edge) > 0)
        {
            printf(" + %.3lf", divided_difference(0, i, deflection_right_edge));
        }
        else if (divided_difference(0, i, deflection_right_edge) < 0)
        {
            printf(" - %.3lf", fabs(divided_difference(0, i, deflection_right_edge)));
        }

        for (int j = 0; j < i; j++)
        {
            printf("(a-%.2lf)", deflection_location[j]);
        }
    }

    puts("");
    puts("");

    printf("Upper Edge: \n");

    printf("f(a) = ");

    for (int i = 0; i < 5; i++)
    {
        if (divided_difference(0, i, deflection_upper_edge) > 0)
        {
            printf(" + %.3lf", divided_difference(0, i, deflection_upper_edge));
        }
        else if (divided_difference(0, i, deflection_upper_edge) < 0)
        {
            printf(" - %.3lf", fabs(divided_difference(0, i, deflection_upper_edge)));
        }

        for (int j = 0; j < i; j++)
        {
            printf("(a-%.2lf)", deflection_location[j]);
        }
    }

    puts("");
    puts("");

    printf("Lower Edge: \n");

    printf("f(a) = ");

    for (int i = 0; i < 5; i++)
    {
        if (divided_difference(0, i, deflection_lower_edge) > 0)
        {
            printf(" + %.3lf", divided_difference(0, i, deflection_lower_edge));
        }
        else if (divided_difference(0, i, deflection_lower_edge) < 0)
        {
            printf(" - %.3lf", fabs(divided_difference(0, i, deflection_lower_edge)));
        }

        for (int j = 0; j < i; j++)
        {
            printf("(a-%.2lf)", deflection_location[j]);
        }
    }

    puts("");
    puts("");
}

void transform_matrix(int original_nodes)
{
    /*
    If the final deflections matrix has n x n values, 
    then we do not know the values for the inner values (the inner (n-2) x (n-2) matrix)
    Thus, we have (n-2)*(n-2) unknowns, constituting a (n-12*(n-2) x (n-2)*(n-2) matrix for Guass-Seidel Method.
    */

    int new_nodes = original_nodes - 2;
    int no_of_unknowns = (int) pow(new_nodes, 2);

    //Initialize Matrix A and B

    A = malloc(no_of_unknowns*sizeof(double *));
    for (int i = 0; i < no_of_unknowns; i++)
    {
        A[i] = malloc(no_of_unknowns*sizeof(double));
    }

    B = malloc(no_of_unknowns*sizeof(double));

    for (int i = 0; i < no_of_unknowns; i++)
    {
        B[i] = 0;

        for (int j = 0; j < no_of_unknowns; j++)
        {
            A[i][j] = 0;
        }
    }

    int row_counter = 0;

    for (int i = 1; i < original_nodes - 1; i++)
    {
        for (int j = 1; j < original_nodes - 1; j++)
        {
            A[row_counter][row_counter] = 4;

            if (final_deflection[j-1][i] != 0) //value to the left
            {
                B[row_counter] += final_deflection[j-1][i];
            }
            else
            {
                A[row_counter][row_counter-1] = -1;
            }

            if (final_deflection[j+1][i] != 0) //value to the right
            {
                B[row_counter] += final_deflection[j+1][i];
            }
            else
            {
                A[row_counter][row_counter+1] = -1;
            }

            if (final_deflection[j][i-1] != 0) //value below
            {
                B[row_counter] += final_deflection[j][i-1];
            }
            else
            {
                A[row_counter][row_counter-new_nodes] = -1;
            }

            if (final_deflection[j][i+1] != 0) //value above
            {
                B[row_counter] += final_deflection[j][i+1];
            }
            else
            {
                A[row_counter][row_counter+new_nodes] = -1;
            }

            row_counter++;
        }
    }

    /*
    puts("");

    for (int i = 0; i < no_of_unknowns; i++)
    {
        for (int j = 0; j < no_of_unknowns; j++)
        {
            printf("%.1lf\t", A[i][j]);
        }
        puts("");
    }

    puts("");

    for (int i = 0; i < no_of_unknowns; i++)
    {
        printf("%.4lf\n", B[i]);
    }

    puts("");
    */
}

double *gauss_seidel(double **A, double *B, int no_of_unknowns, double tolerance)
{
    x = malloc(no_of_unknowns * sizeof(double));

    //Ask for initial guesses
    printf("Enter initial guesses for %d unknown values.\n", no_of_unknowns);

    for (int i = 0; i < no_of_unknowns; i++) 
    {
        printf("x%d: ", i+1);
        scanf("%lf", &x[i]);
    }

    double sum = 0;

    double temp[no_of_unknowns], e[no_of_unknowns];
    for (int i = 0; i < no_of_unknowns; i++)
    {
        e[i] = 1;
    }    

    int checker;

    while (1 == 1) //while true 
    {
        checker = 0;

        for (int i = 0; i < no_of_unknowns; i++)
        {
            for (int j = 0; j < no_of_unknowns; j++)
            {
                if (j != i)
                {
                    sum += A[i][j] * x[j];
                }
            }

            temp[i] = (B[i] - sum) / A[i][i];

            e[i] = fabs(temp[i] - x[i]);
            x[i] = temp[i];
            sum = 0;
        }

        for (int k = 0; k < no_of_unknowns; k++)
        {
            if (e[k] > tolerance)
            {
                checker++;
            }
        }

        //if all errors are less than tolerance, break out of loop
        if (checker == 0)
        {
            break;
        }
    }

    return x;
}

void get_points(int nodes, double spacing, double mew)
{
    x_zero = malloc(sizeof(double));
    y_zero = malloc(sizeof(double));

    no_of_points = 0;

    double x_coor = 0;
    double y_coor = 0;

    //go over rows
    for (int i = 0; i < nodes; i++)
    {
        for (int j = 0; j < nodes - 1; j++)
        {
            if ((final_deflection[j][i] * final_deflection [j+1][i]) < 0) //if two adjacent nodes have opposite signs
            {
                //use linear interpolation
                x_coor = ((spacing / (fabs(final_deflection [j][i]) + fabs(final_deflection[j+1][i]))) * fabs(final_deflection[j][i])) + final_deflection_location[j];
                y_coor = final_deflection_location[i];    

                if ((x_coor > 0.3*mew) && (x_coor < 0.7*mew))
                {
                    no_of_points++;

                    x_zero = realloc(x_zero, no_of_points * sizeof(double));
                    y_zero = realloc(y_zero, no_of_points * sizeof(double));

                    x_zero[no_of_points-1] = x_coor;
                    y_zero[no_of_points-1] = y_coor;
                }                                                   
            }
        }
    }

    //go over columns
    for (int i = 0; i < nodes; i++)
    {
        for (int j = 0; j < nodes - 1; j++)
        {
            if ((final_deflection[i][j] * final_deflection [i][j+1]) < 0) //if two adjacent nodes have opposite signs
            {
                x_coor = final_deflection_location[i];  
                //use linear interpolation
                y_coor = ((spacing / (fabs(final_deflection [i][j]) + fabs(final_deflection[i][j+1]))) * fabs(final_deflection[i][j])) + final_deflection_location[j];

                if ((x_coor > 0.3*mew) && (x_coor < 0.7*mew))
                {
                    no_of_points++;

                    x_zero = realloc(x_zero, no_of_points * sizeof(double));
                    y_zero = realloc(y_zero, no_of_points * sizeof(double));

                    x_zero[no_of_points-1] = x_coor;
                    y_zero[no_of_points-1] = y_coor;
                }
            }
        }
    }
}

double *least_squares(int N, double *x_values, double *y_values)
{
    double *line = malloc(3 * sizeof(double));

    double sum_X = 0;
    for (int i = 0; i < N; i++)
    {
        sum_X += x_values[i];
    }

    double x_mean = sum_X / N;

    double sum_Y = 0;
    for (int i = 0; i < N; i++)
    {
        sum_Y += y_values[i];
    }

    double y_mean = sum_Y / N;

    double sum_XY = 0;
    for (int i = 0; i < N; i++)
    {
        sum_XY += x_values[i] * y_values[i];
    }

    double sum_X2 = 0;
    for (int i = 0; i < N; i++)
    {
        sum_X2 += pow(x_values[i], 2);
    }

    //A
    line[0] = ((N*sum_XY)-(sum_X*sum_Y)) / ((N*sum_X2)-(pow(sum_X,2)));

    //B
    line[1] = y_mean - (line[0]*x_mean);

    double E = 0;
    for (int i = 0; i < N; i++)
    {
        E += pow((line[0]*x_values[i]) + line[1] - y_values[i], 2);
    }

    double E_m = 0;
    for (int i = 0; i < N; i++)
    {
        E_m += pow(y_values[i] - y_mean, 2);
    }

    //r^2
    line[2] = 1 - (E/E_m);

    return line; 
}