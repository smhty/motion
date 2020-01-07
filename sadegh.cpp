
int t_1_a = floor(accel / jerk);
  int t_1_v = floor(sqrt(abs(vel - vInitial) / jerk));
  int t_1_d;
  int t_4;
  int t_2_v, t_2_d;
  double roots[3];
  int num_roots = solveThirdDegree(sign(vel - vInitial) *jerk, sign(vel - vInitial) *(-jerk), 2 * vInitial, -d, roots);
  if (vel >= vInitial)
  {
   t_1_d = floor(max(roots, num_roots));
  }
  else
  {
   if (num_roots <= 1)
   {
    t_1_d = floor(1.5*(1 + sqrt(1 + 6 * vInitial / jerk)));

   }
   else
   {
    t_1_d = floor(roots[1]);
   }
  }
  double tmp[3] = { (double)-t_1_a, (double)-t_1_v, (double)-t_1_d };
  int t_1 = int(-max(tmp, 3));
  int t_2 = 0;
  if (t_1 == 0)
  {
   if (vInitial > 0)
   {
    t_4 = ceil(d / vInitial);
   }
   else
   {
    t_4 = 0;
   }
   j_m = 0;
  }
  else
  {
   if (t_1_a <= MIN(t_1_v, t_1_d))
   {
    t_2_v = floor((abs(vel - vInitial) / (jerk * t_1)) - t_1);
    num_roots = solveThirdDegree(0., sign(vel - vInitial) * jerk * t_1 / 2, sign(vel - vInitial) * 3 * jerk * t_1 * t_1 / 2 - sign(vel - vInitial) * jerk * t_1 + vInitial, sign(vel - vInitial) * jerk * t_1 * t_1 * t_1 - sign(vel - vInitial) * jerk * t_1 * t_1 + 2 * vInitial * t_1 - d, roots);
    if (vel >= vInitial)
    {
     t_2_d = floor(max(roots, num_roots));
    }
    else
    {
     if (num_roots <= 1)
     {
      t_2_d = floor(1 + (vInitial / (jerk * t_1)) - t_1 * 3 / 2);
     }
     else
     {
      t_2_d = floor(roots[0]);
     }
    }
    t_2 = MIN(t_2_v, t_2_d);
   }
   t_4 = ceil((d - (vInitial * (2 * t_1 + t_2) + sign(vel - vInitial) * jerk * t_1 * (t_1 + t_2) * (t_1 - 1 + t_2 / 2))) / (sign(vel - vInitial) * jerk * t_1 * (t_1 + t_2) + vInitial));
   j_m = (d - vInitial * (2 * t_1 + t_2 + t_4)) / (t_1 * (t_1 + t_2) * (t_1 + t_2 / 2 + t_4 - 1));
  }
  a_m = j_m * t_1;
  v_m = a_m * (t_1 + t_2) + vInitial;
  d_m = (v_m - vInitial)* (t_1 + t_2 / 2 + t_4 - 1) + vInitial * (2 * t_1 + t_2 + t_4);
