within ;
package TEGDehydration "Reduced-order dynamic model of a TEG gas dehydration unit"
  extends Modelica.Icons.Package;

  model TEGUnit "TEG dehydration unit: reboiler, surge tank, contactor and rich TEG path"
    extends Modelica.Blocks.Icons.Block;

    // ---- control inputs --------------------------------------------------------
    Modelica.Blocks.Interfaces.RealInput Q_reb(unit="W", start=1.0e5)
      "Reboiler heat duty" annotation (Placement(transformation(extent={{-140,60},{-100,100}})));
    Modelica.Blocks.Interfaces.RealInput m_circ(unit="kg/s", start=0.473)
      "Lean TEG circulation (pump flow, solution basis)"
      annotation (Placement(transformation(extent={{-140,20},{-100,60}})));
    // ---- disturbance inputs ---------------------------------------------------
    Modelica.Blocks.Interfaces.RealInput F_gas(unit="1", start=1.42)
      "Gas flow, million Sm3/d" annotation (Placement(transformation(extent={{-140,-60},{-100,-20}})));
    Modelica.Blocks.Interfaces.RealInput T_gas(unit="1", start=40)
      "Gas temperature at contactor inlet, degC (lean TEG enters 5 K warmer)"
      annotation (Placement(transformation(extent={{-140,-100},{-100,-60}})));

    // ---- outputs ----------------------------------------------------------------
    Modelica.Blocks.Interfaces.RealOutput dryGasWater(unit="1")
      "Dry gas water content, mg/Sm3" annotation (Placement(transformation(extent={{100,70},{120,90}})));
    Modelica.Blocks.Interfaces.RealOutput T_reboiler(unit="1")
      "Reboiler temperature, degC" annotation (Placement(transformation(extent={{100,40},{120,60}})));
    Modelica.Blocks.Interfaces.RealOutput leanTEG(unit="1")
      "Lean TEG purity at the contactor, wt%" annotation (Placement(transformation(extent={{100,10},{120,30}})));
    Modelica.Blocks.Interfaces.RealOutput surgeLevel(unit="1")
      "Surge tank level, 0-1" annotation (Placement(transformation(extent={{100,-20},{120,0}})));
    Modelica.Blocks.Interfaces.RealOutput richTEG(unit="1")
      "Rich TEG purity entering the reboiler, wt%" annotation (Placement(transformation(extent={{100,-50},{120,-30}})));
    Modelica.Blocks.Interfaces.RealOutput boilUp(unit="kg/s")
      "Water vapour leaving the still" annotation (Placement(transformation(extent={{100,-80},{120,-60}})));

    // ---- parameters (equipment and physical properties) -------------------------
    parameter Real cp(unit="J/(kg.K)") = 2600 "Liquid heat capacity, TEG solution";
    parameter Real dHvap(unit="J/kg") = 2.26e6 "Latent heat of water";
    parameter Real T_richIn = 150 "Rich TEG temperature entering the reboiler, degC";
    parameter Real M_weir(unit="kg") = 550 "Reboiler + still holdup at the weir";
    parameter Real k_weir(unit="1/s") = 0.01 "Weir overflow coefficient";
    parameter Real k_evap(unit="kg/(s.K)") = 0.05 "Boil-up per K above bubble point";
    parameter Real M_surgeMax(unit="kg") = 3000 "Surge tank capacity";
    parameter Real tau_contactor(unit="s") = 180 "Contactor liquid mixing time";
    parameter Real tau_rich(unit="s") = 600 "Rich TEG transport: flash drum and exchangers";
    parameter Real dT_smooth = 0.05 "Smoothing width for boil-up and overflow, K or kg";

    // ---- fitted relations (from the calibrated DWSIM model, data/reduced_model_coeffs.json)
    parameter Real cBub[3] = {68.12601092302936, -10.10355851112612, 3.6555717796162437}
      "T_bub(degC) = c1 + c2*ln(w) + c3*ln(w)^2, w water mass fraction, 101.3 kPa";
    parameter Real kMap[10] = {4.127006836105465, -0.8944479518738877, -0.3135933364823048,
      0.26220063000231825, -0.5576886265311166, 0.3662406005743607, -0.002995452135950104,
      -0.4794779481944846, 0.00026420961962964155, 8.892678223368792e-05}
      "ln(dry gas water) polynomial in a = lean wt% - 99, b = ln(ratio/18), c = (T_gas - 40)/5";
    parameter Real dWet[3] = {4.619215685228426, 0.0663578605818007, -0.0001843833886297832}
      "ln(wet gas water, mg/Sm3) = d1 + d2*T_gas + d3*T_gas^2 at 7000 kPa";

    // ---- initialisation -------------------------------------------------------
    parameter Boolean steadyStateInit = true "Start from steady state at the initial inputs";
    parameter Real surgeLevel0 = 0.5 "Initial surge tank level";

    // ---- states ------------------------------------------------------------------
    Real M_reb(unit="kg", start=560) "Reboiler + still liquid mass";
    Real w_reb(start=0.0085, min=1e-4, max=0.2) "Water mass fraction in the reboiler";
    Real T_reb(start=200) "Reboiler temperature, degC";
    Real M_surge(unit="kg", start=1500) "Surge tank liquid mass";
    Real w_surge(start=0.0085, min=1e-4, max=0.2) "Water mass fraction in the surge tank";
    Real w_lean(start=0.0085, min=1e-4, max=0.2) "Water mass fraction of lean TEG in the contactor";
    Real m_rich(unit="kg/s", start=0.49) "Rich TEG flow into the reboiler";
    Real w_rich(start=0.04, min=1e-4, max=0.3) "Water mass fraction of rich TEG into the reboiler";

    // ---- algebraic variables ---------------------------------------------------
    Real T_bub "Bubble temperature of the reboiler liquid, degC";
    Real V(unit="kg/s") "Boil-up (water vapour)";
    Real m_over(unit="kg/s") "Overflow from reboiler to surge tank";
    Real wetGasWater "Wet gas water content, mg/Sm3";
    Real waterIn(unit="kg/s") "Water entering with the gas";
    Real waterAbs(unit="kg/s") "Water absorbed in the contactor";
    Real ratio "Pure TEG / water in wet gas, kg/kg";
    Real a "Lean wt% - 99 (clipped to the map range)";
    Real b "ln(ratio/18) (clipped)";
    Real c "(T_gas - 40)/5";
    Real m_richNow(unit="kg/s") "Rich TEG flow leaving the contactor";
    Real w_richNow "Water fraction of rich TEG leaving the contactor";

  protected
    function softplus "Smooth max(0, x)"
      input Real x;
      input Real eps;
      output Real y;
    algorithm
      y := 0.5*(x + sqrt(x*x + eps*eps));
    end softplus;

  initial equation
    if steadyStateInit then
      der(M_reb) = 0;
      der(w_reb) = 0;
      der(T_reb) = 0;
      der(w_surge) = 0;
      der(w_lean) = 0;
      der(m_rich) = 0;
      der(w_rich) = 0;
    else
      M_reb = M_weir + 10;
      w_reb = 0.0085;
      T_reb = 200;
      w_surge = 0.0085;
      w_lean = 0.0085;
      m_rich = 0.49;
      w_rich = 0.04;
    end if;
    M_surge = surgeLevel0*M_surgeMax;

  equation
    // reboiler + still: liquid boils off water at its bubble point
    T_bub = cBub[1] + cBub[2]*log(w_reb) + cBub[3]*log(w_reb)^2;
    V = k_evap*softplus(T_reb - T_bub, dT_smooth);
    m_over = k_weir*softplus(M_reb - M_weir, dT_smooth);
    der(M_reb) = m_rich - m_over - V;
    M_reb*der(w_reb) = m_rich*(w_rich - w_reb) - V*(1 - w_reb);
    M_reb*cp*der(T_reb) = m_rich*cp*(T_richIn - T_reb) + Q_reb - V*dHvap;

    // surge tank (well mixed); the pump draws m_circ
    der(M_surge) = m_over - m_circ;
    M_surge*der(w_surge) = m_over*(w_reb - w_surge);

    // contactor: liquid holdup mixing, then steady performance map
    tau_contactor*der(w_lean) = w_surge - w_lean;
    wetGasWater = exp(dWet[1] + dWet[2]*T_gas + dWet[3]*T_gas^2);
    waterIn = F_gas*1e6/86400*wetGasWater*1e-6;
    ratio = m_circ*(1 - w_lean)/max(waterIn, 1e-6);
    a = min(max(100*(1 - w_lean) - 99, -1.0), 0.9);
    b = log(min(max(ratio, 6), 45)/18);
    c = (T_gas - 40)/5;
    dryGasWater = min(exp(kMap[1] + kMap[2]*a + kMap[3]*b + kMap[4]*c + kMap[5]*a^2 + kMap[6]*b^2
      + kMap[7]*c^2 + kMap[8]*a*b + kMap[9]*a*c + kMap[10]*b*c), wetGasWater);
    waterAbs = F_gas*1e6/86400*(wetGasWater - dryGasWater)*1e-6;
    m_richNow = m_circ + waterAbs;
    w_richNow = (m_circ*w_lean + waterAbs)/m_richNow;

    // rich TEG transport lag (flash drum, lean/rich exchanger)
    tau_rich*der(m_rich) = m_richNow - m_rich;
    tau_rich*der(w_rich) = w_richNow - w_rich;

    // outputs
    T_reboiler = T_reb;
    leanTEG = 100*(1 - w_lean);
    surgeLevel = M_surge/M_surgeMax;
    richTEG = 100*(1 - w_rich);
    boilUp = V;
    annotation (Documentation(info="<html>
<p>Reduced-order dynamic model of a TEG dehydration unit for a 1.42 million Sm3/d coal-seam gas stream at
7000 kPa. Dynamic: reboiler/still holdup, composition and energy; surge tank level and mixing; contactor liquid
mixing; rich TEG transport lag. Steady: contactor performance map and TEG-water bubble curve, both fitted to a
DWSIM model with Peng-Robinson parameters calibrated on public data (see the repository README).</p>
<p>Inputs: Q_reb (W), m_circ (kg/s); disturbances F_gas (million Sm3/d), T_gas (degC).
Outputs: dry gas water (mg/Sm3, spec 112), reboiler temperature, lean and rich TEG purity, surge level, boil-up.</p>
<p>Map validity: lean TEG 98.5-99.7 wt%, TEG/water ratio 8-36 kg/kg, gas 30-45 degC. Equipment holdups and time
constants are assumed design values.</p>
</html>"));
  end TEGUnit;

  package Examples "Test cases"
    extends Modelica.Icons.ExamplesPackage;

    model StepTest "Open-loop steps: reboiler duty, circulation, gas temperature, gas flow"
      extends Modelica.Icons.Example;
      TEGUnit unit annotation (Placement(transformation(extent={{-10,-10},{10,10}})));
      Modelica.Blocks.Sources.Step Q(height=-1.5e4, offset=1.0e5, startTime=2*3600)
        "Reboiler duty: -15 kW at 2 h" annotation (Placement(transformation(extent={{-80,60},{-60,80}})));
      Modelica.Blocks.Sources.Step circ(height=-0.15, offset=0.473, startTime=8*3600)
        "Circulation: -0.15 kg/s at 8 h" annotation (Placement(transformation(extent={{-80,20},{-60,40}})));
      Modelica.Blocks.Sources.Step gasFlow(height=0.3, offset=1.42, startTime=14*3600)
        "Gas flow: +0.3 million Sm3/d at 14 h" annotation (Placement(transformation(extent={{-80,-40},{-60,-20}})));
      Modelica.Blocks.Sources.Step temp(height=5, offset=40, startTime=20*3600)
        "Gas temperature: +5 K at 20 h" annotation (Placement(transformation(extent={{-80,-80},{-60,-60}})));
    equation
      connect(Q.y, unit.Q_reb) annotation (Line(points={{-59,70},{-30,70},{-30,8},{-12,8}}, color={0,0,127}));
      connect(circ.y, unit.m_circ) annotation (Line(points={{-59,30},{-40,30},{-40,4},{-12,4}}, color={0,0,127}));
      connect(gasFlow.y, unit.F_gas) annotation (Line(points={{-59,-30},{-40,-30},{-40,-4},{-12,-4}}, color={0,0,127}));
      connect(temp.y, unit.T_gas) annotation (Line(points={{-59,-70},{-30,-70},{-30,-8},{-12,-8}}, color={0,0,127}));
      annotation (experiment(StopTime=93600, Interval=60, Tolerance=1e-06));
    end StepTest;
  end Examples;
  annotation (uses(Modelica(version="4.0.0")));
end TEGDehydration;
