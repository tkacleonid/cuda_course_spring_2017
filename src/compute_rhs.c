#include <math.h>
#include "header.h"

void compute_rhs()
{
    int i, j, k, m;
    double aux, rho_inv, uijk, up1, um1, vijk, vp1, vm1, wijk, wp1, wm1;


    if (timeron) timer_start(t_rhs);
    for (k = 0; k <= nz - 1; k++) 
    {
        for (j = 0; j <= ny - 1; j++) 
        {
            for (i = 0; i <= nx - 1; i++) 
            {
                rho_inv = 1.0 / u[k][j][i][0];
                rho_i[k][j][i] = rho_inv;
                us[k][j][i] = u[k][j][i][1] * rho_inv;
                vs[k][j][i] = u[k][j][i][2] * rho_inv;
                ws[k][j][i] = u[k][j][i][3] * rho_inv;
                square[k][j][i] = 0.5* (u[k][j][i][1] * u[k][j][i][1] + u[k][j][i][2] * u[k][j][i][2] + u[k][j][i][3] * u[k][j][i][3]) * rho_inv;
                qs[k][j][i] = square[k][j][i] * rho_inv;
                aux = c1c2*rho_inv* (u[k][j][i][4] - square[k][j][i]);
                speed[k][j][i] = sqrt(aux);

                for (m = 0; m < 5; m++)
                    rhs[k][j][i][m] = forcing[k][j][i][m];
            }
        }
    }
    

    if (timeron) timer_start(t_rhsx);
    for (k = 1; k <= nz2; k++) 
    {
        for (j = 1; j <= ny2; j++) 
        {
            for (i = 1; i <= nx2; i++) 
            {
                uijk = us[k][j][i];
                up1 = us[k][j][i + 1];
                um1 = us[k][j][i - 1];

                rhs[k][j][i][0] = rhs[k][j][i][0] + dx1tx1 *
                    (u[k][j][i + 1][0] - 2.0*u[k][j][i][0] + u[k][j][i - 1][0]) -
                    tx2 * (u[k][j][i + 1][1] - u[k][j][i - 1][1]);

                rhs[k][j][i][1] = rhs[k][j][i][1] + dx2tx1 *
                    (u[k][j][i + 1][1] - 2.0*u[k][j][i][1] + u[k][j][i - 1][1]) +
                    xxcon2*con43 * (up1 - 2.0*uijk + um1) -
                    tx2 * (u[k][j][i + 1][1] * up1 - u[k][j][i - 1][1] * um1 +
                    (u[k][j][i + 1][4] - square[k][j][i + 1] -
                    u[k][j][i - 1][4] + square[k][j][i - 1]) * c2);

                rhs[k][j][i][2] = rhs[k][j][i][2] + dx3tx1 *
                    (u[k][j][i + 1][2] - 2.0*u[k][j][i][2] + u[k][j][i - 1][2]) +
                    xxcon2 * (vs[k][j][i + 1] - 2.0*vs[k][j][i] + vs[k][j][i - 1]) -
                    tx2 * (u[k][j][i + 1][2] * up1 - u[k][j][i - 1][2] * um1);

                rhs[k][j][i][3] = rhs[k][j][i][3] + dx4tx1 *
                    (u[k][j][i + 1][3] - 2.0*u[k][j][i][3] + u[k][j][i - 1][3]) +
                    xxcon2 * (ws[k][j][i + 1] - 2.0*ws[k][j][i] + ws[k][j][i - 1]) -
                    tx2 * (u[k][j][i + 1][3] * up1 - u[k][j][i - 1][3] * um1);

                rhs[k][j][i][4] = rhs[k][j][i][4] + dx5tx1 *
                    (u[k][j][i + 1][4] - 2.0*u[k][j][i][4] + u[k][j][i - 1][4]) +
                    xxcon3 * (qs[k][j][i + 1] - 2.0*qs[k][j][i] + qs[k][j][i - 1]) +
                    xxcon4 * (up1*up1 - 2.0*uijk*uijk + um1*um1) +
                    xxcon5 * (u[k][j][i + 1][4] * rho_i[k][j][i + 1] -
                    2.0*u[k][j][i][4] * rho_i[k][j][i] +
                    u[k][j][i - 1][4] * rho_i[k][j][i - 1]) -
                    tx2 * ((c1*u[k][j][i + 1][4] - c2*square[k][j][i + 1])*up1 -
                    (c1*u[k][j][i - 1][4] - c2*square[k][j][i - 1])*um1);

                if (i == 1)
                {
                    for (m = 0; m < 5; m++)
                        rhs[k][j][i][m] = rhs[k][j][i][m] - dssp * (5.0*u[k][j][i][m] - 4.0*u[k][j][i + 1][m] + u[k][j][i + 2][m]);
                }
                else if (i == 2)
                {
                    for (m = 0; m < 5; m++)
                        rhs[k][j][i][m] = rhs[k][j][i][m] - dssp * (-4.0*u[k][j][i - 1][m] + 6.0*u[k][j][i][m] - 4.0*u[k][j][i + 1][m] + u[k][j][i + 2][m]);
                }
                else if (i == nx2 - 1)
                {
                    for (m = 0; m < 5; m++)
                        rhs[k][j][i][m] = rhs[k][j][i][m] - dssp * (u[k][j][i - 2][m] - 4.0*u[k][j][i - 1][m] + 6.0*u[k][j][i][m] - 4.0*u[k][j][i + 1][m]);
                }
                else if (i == nx2)
                {
                    for (m = 0; m < 5; m++)
                        rhs[k][j][i][m] = rhs[k][j][i][m] - dssp * (u[k][j][i - 2][m] - 4.0*u[k][j][i - 1][m] + 5.0*u[k][j][i][m]);
                }
                else
                {
                    for (m = 0; m < 5; m++)
                        rhs[k][j][i][m] = rhs[k][j][i][m] - dssp * (u[k][j][i - 2][m] - 4.0*u[k][j][i - 1][m] + 6.0*u[k][j][i][m] - 4.0*u[k][j][i + 1][m] + u[k][j][i + 2][m]);
                }
            }
        }
    }
    if (timeron) timer_stop(t_rhsx);


    if (timeron) timer_start(t_rhsy);
    for (k = 1; k <= nz2; k++) 
    {
        for (j = 1; j <= ny2; j++) 
        {
            for (i = 1; i <= nx2; i++) 
            {
                vijk = vs[k][j][i];
                vp1 = vs[k][j + 1][i];
                vm1 = vs[k][j - 1][i];

                rhs[k][j][i][0] = rhs[k][j][i][0] + dy1ty1 *
                    (u[k][j + 1][i][0] - 2.0*u[k][j][i][0] + u[k][j - 1][i][0]) -
                    ty2 * (u[k][j + 1][i][2] - u[k][j - 1][i][2]);

                rhs[k][j][i][1] = rhs[k][j][i][1] + dy2ty1 *
                    (u[k][j + 1][i][1] - 2.0*u[k][j][i][1] + u[k][j - 1][i][1]) +
                    yycon2 * (us[k][j + 1][i] - 2.0*us[k][j][i] + us[k][j - 1][i]) -
                    ty2 * (u[k][j + 1][i][1] * vp1 - u[k][j - 1][i][1] * vm1);

                rhs[k][j][i][2] = rhs[k][j][i][2] + dy3ty1 *
                    (u[k][j + 1][i][2] - 2.0*u[k][j][i][2] + u[k][j - 1][i][2]) +
                    yycon2*con43 * (vp1 - 2.0*vijk + vm1) -
                    ty2 * (u[k][j + 1][i][2] * vp1 - u[k][j - 1][i][2] * vm1 +
                    (u[k][j + 1][i][4] - square[k][j + 1][i] -
                    u[k][j - 1][i][4] + square[k][j - 1][i]) * c2);

                rhs[k][j][i][3] = rhs[k][j][i][3] + dy4ty1 *
                    (u[k][j + 1][i][3] - 2.0*u[k][j][i][3] + u[k][j - 1][i][3]) +
                    yycon2 * (ws[k][j + 1][i] - 2.0*ws[k][j][i] + ws[k][j - 1][i]) -
                    ty2 * (u[k][j + 1][i][3] * vp1 - u[k][j - 1][i][3] * vm1);

                rhs[k][j][i][4] = rhs[k][j][i][4] + dy5ty1 *
                    (u[k][j + 1][i][4] - 2.0*u[k][j][i][4] + u[k][j - 1][i][4]) +
                    yycon3 * (qs[k][j + 1][i] - 2.0*qs[k][j][i] + qs[k][j - 1][i]) +
                    yycon4 * (vp1*vp1 - 2.0*vijk*vijk + vm1*vm1) +
                    yycon5 * (u[k][j + 1][i][4] * rho_i[k][j + 1][i] -
                    2.0*u[k][j][i][4] * rho_i[k][j][i] +
                    u[k][j - 1][i][4] * rho_i[k][j - 1][i]) -
                    ty2 * ((c1*u[k][j + 1][i][4] - c2*square[k][j + 1][i]) * vp1 -
                    (c1*u[k][j - 1][i][4] - c2*square[k][j - 1][i]) * vm1);

                if (j == 1)
                {
                    for (m = 0; m < 5; m++)
                        rhs[k][j][i][m] = rhs[k][j][i][m] - dssp * (5.0*u[k][j][i][m] - 4.0*u[k][j + 1][i][m] + u[k][j + 2][i][m]);
                }
                else if (j == 2)
                {
                    for (m = 0; m < 5; m++)
                        rhs[k][j][i][m] = rhs[k][j][i][m] - dssp * (-4.0*u[k][j - 1][i][m] + 6.0*u[k][j][i][m] - 4.0*u[k][j + 1][i][m] + u[k][j + 2][i][m]);
                }
                else if (j == ny2 - 1)
                {
                    for (m = 0; m < 5; m++)
                        rhs[k][j][i][m] = rhs[k][j][i][m] - dssp * (u[k][j - 2][i][m] - 4.0*u[k][j - 1][i][m] + 6.0*u[k][j][i][m] - 4.0*u[k][j + 1][i][m]);
                }
                else if (j == ny2)
                {
                    for (m = 0; m < 5; m++)
                        rhs[k][j][i][m] = rhs[k][j][i][m] - dssp * (u[k][j - 2][i][m] - 4.0*u[k][j - 1][i][m] + 5.0*u[k][j][i][m]);
                }
                else
                {
                    for (m = 0; m < 5; m++)
                        rhs[k][j][i][m] = rhs[k][j][i][m] - dssp * (u[k][j - 2][i][m] - 4.0*u[k][j - 1][i][m] + 6.0*u[k][j][i][m] - 4.0*u[k][j + 1][i][m] + u[k][j + 2][i][m]);
                }
            }
        }
    }
    if (timeron) timer_stop(t_rhsy);


    if (timeron) timer_start(t_rhsz);
    for (k = 1; k <= nz2; k++) 
    {
        for (j = 1; j <= ny2; j++) 
        {
            for (i = 1; i <= nx2; i++)
            {
                wijk = ws[k][j][i];
                wp1 = ws[k + 1][j][i];
                wm1 = ws[k - 1][j][i];

                rhs[k][j][i][0] = rhs[k][j][i][0] + dz1tz1 *
                    (u[k + 1][j][i][0] - 2.0*u[k][j][i][0] + u[k - 1][j][i][0]) -
                    tz2 * (u[k + 1][j][i][3] - u[k - 1][j][i][3]);

                rhs[k][j][i][1] = rhs[k][j][i][1] + dz2tz1 *
                    (u[k + 1][j][i][1] - 2.0*u[k][j][i][1] + u[k - 1][j][i][1]) +
                    zzcon2 * (us[k + 1][j][i] - 2.0*us[k][j][i] + us[k - 1][j][i]) -
                    tz2 * (u[k + 1][j][i][1] * wp1 - u[k - 1][j][i][1] * wm1);

                rhs[k][j][i][2] = rhs[k][j][i][2] + dz3tz1 *
                    (u[k + 1][j][i][2] - 2.0*u[k][j][i][2] + u[k - 1][j][i][2]) +
                    zzcon2 * (vs[k + 1][j][i] - 2.0*vs[k][j][i] + vs[k - 1][j][i]) -
                    tz2 * (u[k + 1][j][i][2] * wp1 - u[k - 1][j][i][2] * wm1);

                rhs[k][j][i][3] = rhs[k][j][i][3] + dz4tz1 *
                    (u[k + 1][j][i][3] - 2.0*u[k][j][i][3] + u[k - 1][j][i][3]) +
                    zzcon2*con43 * (wp1 - 2.0*wijk + wm1) -
                    tz2 * (u[k + 1][j][i][3] * wp1 - u[k - 1][j][i][3] * wm1 +
                    (u[k + 1][j][i][4] - square[k + 1][j][i] -
                    u[k - 1][j][i][4] + square[k - 1][j][i]) * c2);

                rhs[k][j][i][4] = rhs[k][j][i][4] + dz5tz1 *
                    (u[k + 1][j][i][4] - 2.0*u[k][j][i][4] + u[k - 1][j][i][4]) +
                    zzcon3 * (qs[k + 1][j][i] - 2.0*qs[k][j][i] + qs[k - 1][j][i]) +
                    zzcon4 * (wp1*wp1 - 2.0*wijk*wijk + wm1*wm1) +
                    zzcon5 * (u[k + 1][j][i][4] * rho_i[k + 1][j][i] -
                    2.0*u[k][j][i][4] * rho_i[k][j][i] +
                    u[k - 1][j][i][4] * rho_i[k - 1][j][i]) -
                    tz2 * ((c1*u[k + 1][j][i][4] - c2*square[k + 1][j][i])*wp1 -
                    (c1*u[k - 1][j][i][4] - c2*square[k - 1][j][i])*wm1);

                if (k == 1)
                {
                    for (m = 0; m < 5; m++)
                        rhs[k][j][i][m] = rhs[k][j][i][m] - dssp * (5.0*u[k][j][i][m] - 4.0*u[k + 1][j][i][m] + u[k + 2][j][i][m]);
                }
                else if (k == 2)
                {
                    for (m = 0; m < 5; m++)
                        rhs[k][j][i][m] = rhs[k][j][i][m] - dssp * (-4.0*u[k - 1][j][i][m] + 6.0*u[k][j][i][m] - 4.0*u[k + 1][j][i][m] + u[k + 2][j][i][m]);
                }
                else if (k == nz2 - 1)
                {
                    for (m = 0; m < 5; m++)
                        rhs[k][j][i][m] = rhs[k][j][i][m] - dssp * (u[k - 2][j][i][m] - 4.0*u[k - 1][j][i][m] + 6.0*u[k][j][i][m] - 4.0*u[k + 1][j][i][m]);
                }
                else if (k == nz2)
                {
                    for (m = 0; m < 5; m++)
                        rhs[k][j][i][m] = rhs[k][j][i][m] - dssp * (u[k - 2][j][i][m] - 4.0*u[k - 1][j][i][m] + 5.0*u[k][j][i][m]);
                }
                else
                {
                    for (m = 0; m < 5; m++)
                        rhs[k][j][i][m] = rhs[k][j][i][m] - dssp * (u[k - 2][j][i][m] - 4.0*u[k - 1][j][i][m] + 6.0*u[k][j][i][m] - 4.0*u[k + 1][j][i][m] + u[k + 2][j][i][m]);
                }
            }
        }
    }   
    if (timeron) timer_stop(t_rhsz);

    for (k = 1; k <= nz2; k++) 
    {
        for (j = 1; j <= ny2; j++) 
        {
            for (i = 1; i <= nx2; i++) 
            {
                for (m = 0; m < 5; m++)
                    rhs[k][j][i][m] = rhs[k][j][i][m] * dt;
            }
        }
    }
    if (timeron) timer_stop(t_rhs);
}
