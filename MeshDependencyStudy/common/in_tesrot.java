// Copyright 2020 Sebastian Engel
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
// THE SOFTWARE.
//
//
//
//
//
// #############################################################################
// Java-Macro for Siemens StarCCM+
//
// Part of the Work:
// Multi-objective Design Optimization of a Tesla-type Rotary Blood Pump
//
// sets up the simulation and runs it
//
// written by Sebastian Engel
// #############################################################################

import java.util.*;

import star.turbulence.*;
import star.kwturb.*;
import star.material.*;
import star.common.*;
import star.base.neo.*;
import star.vis.*;
import star.base.report.*;
import star.flow.*;
import star.coupledflow.*;
import star.passivescalar.*;
import star.segregatedflow.*;
import star.metrics.*;
import star.energy.*;
import star.walldistance.*;

public class in_tesrot extends StarMacro {

    public void execute() {

        prep();
        solve();
        save();

    }

    private void prep() {

        Simulation sim = getActiveSimulation();

        // Importing Mesh
        ImportManager importManager = sim.getImportManager();
        importManager.importMeshFiles(new StringVector(new String[]{resolvePath("slice.msh")}), NeoProperty.fromString("{\'FileOptions\': [{\'Sequence\': 51}]}"));
        Region region_fluid = sim.getRegionManager().getRegion("FLUID");

        // Setting Physics
        PhysicsContinuum physicsContinuum_fluid = ((PhysicsContinuum) sim.getContinuumManager().getContinuum("Physics 1"));

        physicsContinuum_fluid.enable(SingleComponentLiquidModel.class);

        // physicsContinuum_fluid.enable(CoupledFlowModel.class);
        // CoupledImplicitSolver coupledImplicitSolver_0 = ((CoupledImplicitSolver) sim.getSolverManager().getSolver(CoupledImplicitSolver.class));
        // coupledImplicitSolver_0.setLeaveTemporaryStorage(true);
        // CoupledFlowModel coupledFlowModel_0 = physicsContinuum_fluid.getModelManager().getModel(CoupledFlowModel.class);
        // coupledFlowModel_0.getMaxUref().setValue(100.0);		
        physicsContinuum_fluid.enable(SegregatedFlowModel.class);
        SegregatedFlowSolver segregatedFlowSolver_0 = ((SegregatedFlowSolver) sim.getSolverManager().getSolver(SegregatedFlowSolver.class));
        segregatedFlowSolver_0.setLeaveTemporaryStorage(true);

        physicsContinuum_fluid.enable(SteadyModel.class);

        physicsContinuum_fluid.enable(ConstantDensityModel.class);
        physicsContinuum_fluid.enable(TurbulentModel.class);
        physicsContinuum_fluid.enable(RansTurbulenceModel.class);
        physicsContinuum_fluid.enable(KOmegaTurbulence.class);
        physicsContinuum_fluid.enable(SstKwTurbModel.class);
        physicsContinuum_fluid.enable(ThreeDimensionalModel.class);

        physicsContinuum_fluid.enable(KwAllYplusWallTreatment.class);
        physicsContinuum_fluid.enable(GammaReThetaTransitionModel.class);
        KwAllYplusWallTreatment kwAllYplusWallTreatment_0 = physicsContinuum_fluid.getModelManager().getModel(KwAllYplusWallTreatment.class);
        physicsContinuum_fluid.disableModel(kwAllYplusWallTreatment_0);
        physicsContinuum_fluid.enable(KwLowYplusWallTreatment.class);

        SingleComponentLiquidModel singleComponentLiquidModel = physicsContinuum_fluid.getModelManager().getModel(SingleComponentLiquidModel.class);
        Liquid liquid_blood = ((Liquid) singleComponentLiquidModel.getMaterial());
        liquid_blood.setPresentationName("Blood");

        ConstantMaterialPropertyMethod constantMaterialPropertyMethod_bloodDensity = ((ConstantMaterialPropertyMethod) liquid_blood.getMaterialProperties().getMaterialProperty(ConstantDensityProperty.class).getMethod());
        constantMaterialPropertyMethod_bloodDensity.getQuantity().setValue(1059.0);

        liquid_blood.getMaterialProperties().getMaterialProperty(DynamicViscosityProperty.class).setMethod(NonNewtonianGeneralizedCarreauYasudaMethod.class);

        NonNewtonianGeneralizedCarreauYasudaMethod CarreauYasudaMethod_blood = ((NonNewtonianGeneralizedCarreauYasudaMethod) liquid_blood.getMaterialProperties().getMaterialProperty(DynamicViscosityProperty.class).getMethod());
        CarreauYasudaMethod_blood.setPowerConstant(0.392);
        CarreauYasudaMethod_blood.setAParameter(0.644);
        CarreauYasudaMethod_blood.getZeroShearViscosity().setValue(0.022);
        CarreauYasudaMethod_blood.getInfiniteShearViscosity().setValue(0.0022);
        CarreauYasudaMethod_blood.getRelaxationTime().setValue(0.11);

        // Passive Scalar settings for Hemolysis Ratio
        physicsContinuum_fluid.enable(PassiveScalarModel.class);
        PassiveScalarModel passiveScalarModel = physicsContinuum_fluid.getModelManager().getModel(PassiveScalarModel.class);
        PassiveScalarMaterial passiveScalarMaterial_hemolysisRatio = passiveScalarModel.getPassiveScalarManager().createPassiveScalarMaterial(PassiveScalarMaterial.class);
        passiveScalarMaterial_hemolysisRatio.setPresentationName("HemolysisRatio");
        passiveScalarMaterial_hemolysisRatio.getTransportOption().setSelected(PassiveScalarTransportOption.CONVECTION_ONLY);
        passiveScalarMaterial_hemolysisRatio.setMaxAllowable(1.0);
        passiveScalarMaterial_hemolysisRatio.getClipMode().setSelected(PassiveScalarClipMode.CLIP_BOTH);

        Boundary boundary_per0 = region_fluid.getBoundaryManager().getBoundary("SYM2");
        Boundary boundary_per1 = region_fluid.getBoundaryManager().getBoundary("SYM1");

        // BoundaryInterface boundaryInterface_periodic = sim.getInterfaceManager().createBoundaryInterface(boundary_per0, boundary_per1, "Interface", true);
        // boundaryInterface_periodic.getTopology().setSelected(InterfaceConfigurationOption.PERIODIC);
        // InterfacePeriodicTransformSpecification interfacePeriodicTransformSpecification_0 = boundaryInterface_periodic.getPeriodicTransform();
        // interfacePeriodicTransformSpecification_0.getPeriodicityOption().setSelected(PeriodicityOption.TRANSLATION);
        // InterfaceToleranceCondition interfaceToleranceCondition_0 = boundaryInterface_periodic.getValues().get(InterfaceToleranceCondition.class);
        // interfaceToleranceCondition_0.setTolerance(5.0E-6);		
        Boundary boundary_disks = region_fluid.getBoundaryManager().getBoundary("DISKS");
        boundary_disks.getConditions().get(WallSlidingOption.class).setSelected(WallSlidingOption.ROTATION_RATE);
        WallRelativeRotationProfile wallRelativeRotationProfile_disks = boundary_disks.getValues().get(WallRelativeRotationProfile.class);
        wallRelativeRotationProfile_disks.setMethod(FunctionScalarProfileMethod.class);

        wallRelativeRotationProfile_disks.getMethod(ConstantScalarProfileMethod.class).getQuantity().setValue(300.0);

        IntermittencyProfile intermittencyProfile_0 = physicsContinuum_fluid.getInitialConditions().get(IntermittencyProfile.class);
        intermittencyProfile_0.getMethod(ConstantScalarProfileMethod.class).getQuantity().setValue(0.0);

        PassiveScalarProfile passiveScalarProfile_0 = physicsContinuum_fluid.getInitialConditions().get(PassiveScalarProfile.class);
        passiveScalarProfile_0.getMethod(ConstantArrayProfileMethod.class).getQuantity().setArray(new DoubleVector(new double[]{1.0E-15}));

        InitialPressureProfile initialPressureProfile_0 = physicsContinuum_fluid.getInitialConditions().get(InitialPressureProfile.class);
        initialPressureProfile_0.getMethod(ConstantScalarProfileMethod.class).getQuantity().setValue(-14000.0);

        TurbulentViscosityRatioProfile turbulentViscosityRatioProfile_0 = physicsContinuum_fluid.getInitialConditions().get(TurbulentViscosityRatioProfile.class);
        turbulentViscosityRatioProfile_0.getMethod(ConstantScalarProfileMethod.class).getQuantity().setValue(3.0);

        VelocitySolver velocitySolver_0 = segregatedFlowSolver_0.getVelocitySolver();
        velocitySolver_0.setUrf(0.8);

        PressureSolver pressureSolver_0 = segregatedFlowSolver_0.getPressureSolver();
        pressureSolver_0.setUrf(0.2);

        KwTurbSolver kwTurbSolver_0 = ((KwTurbSolver) sim.getSolverManager().getSolver(KwTurbSolver.class));
        kwTurbSolver_0.setLeaveTemporaryStorage(true);
        kwTurbSolver_0.setUrf(1.0);

        GammaReThetaTransitionSolver gammaReThetaTransitionSolver_0 = ((GammaReThetaTransitionSolver) sim.getSolverManager().getSolver(GammaReThetaTransitionSolver.class));
        gammaReThetaTransitionSolver_0.setUrf(1.0);

        PassiveScalarSolver passiveScalarSolver_0 = ((PassiveScalarSolver) sim.getSolverManager().getSolver(PassiveScalarSolver.class));
        passiveScalarSolver_0.setUrf(1.0);

        // Stopping Criteria
        ResidualMonitor residualMonitor_xMomentum = ((ResidualMonitor) sim.getMonitorManager().getMonitor("X-momentum"));
        ResidualMonitor residualMonitor_yMomentum = ((ResidualMonitor) sim.getMonitorManager().getMonitor("Y-momentum"));
        ResidualMonitor residualMonitor_zMomentum = ((ResidualMonitor) sim.getMonitorManager().getMonitor("Z-momentum"));
        ResidualMonitor residualMonitor_continuity = ((ResidualMonitor) sim.getMonitorManager().getMonitor("Continuity"));
        ResidualMonitor residualMonitor_tke = ((ResidualMonitor) sim.getMonitorManager().getMonitor("Tke"));
        ResidualMonitor residualMonitor_sdr = ((ResidualMonitor) sim.getMonitorManager().getMonitor("Sdr"));
        ResidualMonitor residualMonitor_rtheta = ((ResidualMonitor) sim.getMonitorManager().getMonitor("ReTheta_t"));
        ResidualMonitor residualMonitor_intermitt = ((ResidualMonitor) sim.getMonitorManager().getMonitor("Intermittency"));

        MonitorIterationStoppingCriterion stopCrit_xMomentum = residualMonitor_xMomentum.createIterationStoppingCriterion();
        MonitorIterationStoppingCriterion stopCrit_yMomentum = residualMonitor_yMomentum.createIterationStoppingCriterion();
        MonitorIterationStoppingCriterion stopCrit_zMomentum = residualMonitor_zMomentum.createIterationStoppingCriterion();
        MonitorIterationStoppingCriterion stopCrit_continuity = residualMonitor_continuity.createIterationStoppingCriterion();
        MonitorIterationStoppingCriterion stopCrit_tke = residualMonitor_tke.createIterationStoppingCriterion();
        MonitorIterationStoppingCriterion stopCrit_sdr = residualMonitor_sdr.createIterationStoppingCriterion();
        MonitorIterationStoppingCriterion stopCrit_rtheta = residualMonitor_rtheta.createIterationStoppingCriterion();
        MonitorIterationStoppingCriterion stopCrit_intermitt = residualMonitor_intermitt.createIterationStoppingCriterion();

        stopCrit_xMomentum.getLogicalOption().setSelected(SolverStoppingCriterionLogicalOption.AND);
        stopCrit_yMomentum.getLogicalOption().setSelected(SolverStoppingCriterionLogicalOption.AND);
        stopCrit_zMomentum.getLogicalOption().setSelected(SolverStoppingCriterionLogicalOption.AND);
        stopCrit_continuity.getLogicalOption().setSelected(SolverStoppingCriterionLogicalOption.AND);
        stopCrit_tke.getLogicalOption().setSelected(SolverStoppingCriterionLogicalOption.AND);
        stopCrit_sdr.getLogicalOption().setSelected(SolverStoppingCriterionLogicalOption.AND);
        stopCrit_rtheta.getLogicalOption().setSelected(SolverStoppingCriterionLogicalOption.AND);
        stopCrit_intermitt.getLogicalOption().setSelected(SolverStoppingCriterionLogicalOption.AND);

        MonitorIterationStoppingCriterionMinLimitType stopCritLimit_0 = ((MonitorIterationStoppingCriterionMinLimitType) stopCrit_xMomentum.getCriterionType());
        MonitorIterationStoppingCriterionMinLimitType stopCritLimit_1 = ((MonitorIterationStoppingCriterionMinLimitType) stopCrit_yMomentum.getCriterionType());
        MonitorIterationStoppingCriterionMinLimitType stopCritLimit_2 = ((MonitorIterationStoppingCriterionMinLimitType) stopCrit_zMomentum.getCriterionType());
        MonitorIterationStoppingCriterionMinLimitType stopCritLimit_3 = ((MonitorIterationStoppingCriterionMinLimitType) stopCrit_continuity.getCriterionType());
        MonitorIterationStoppingCriterionMinLimitType stopCritLimit_4 = ((MonitorIterationStoppingCriterionMinLimitType) stopCrit_tke.getCriterionType());
        MonitorIterationStoppingCriterionMinLimitType stopCritLimit_5 = ((MonitorIterationStoppingCriterionMinLimitType) stopCrit_sdr.getCriterionType());
        MonitorIterationStoppingCriterionMinLimitType stopCritLimit_6 = ((MonitorIterationStoppingCriterionMinLimitType) stopCrit_rtheta.getCriterionType());
        MonitorIterationStoppingCriterionMinLimitType stopCritLimit_7 = ((MonitorIterationStoppingCriterionMinLimitType) stopCrit_intermitt.getCriterionType());

        StepStoppingCriterion maxStepStoppingCriterion
                = ((StepStoppingCriterion) sim.getSolverStoppingCriterionManager().getSolverStoppingCriterion("Maximum Steps"));
        maxStepStoppingCriterion.setIsUsed(true);
        maxStepStoppingCriterion.setMaximumNumberSteps(10000);

        stopCritLimit_0.getLimit().setValue(1.0E-3);
        stopCritLimit_1.getLimit().setValue(1.0E-3);
        stopCritLimit_2.getLimit().setValue(1.0E-3);
        stopCritLimit_3.getLimit().setValue(1.0E-3);
        stopCritLimit_4.getLimit().setValue(1.0E-3);
        stopCritLimit_5.getLimit().setValue(1.0E-3);
        stopCritLimit_6.getLimit().setValue(1.0E-3);
        stopCritLimit_7.getLimit().setValue(1.0E-3);

        // Setting up Field Functions
        UserFieldFunction uFF_GWparam = sim.getFieldFunctionManager().createFieldFunction();
        UserFieldFunction uFF_massFlow = sim.getFieldFunctionManager().createFieldFunction();
        UserFieldFunction uFF_thresholdValue = sim.getFieldFunctionManager().createFieldFunction();
        UserFieldFunction uFF_equivStress = sim.getFieldFunctionManager().createFieldFunction();
        UserFieldFunction uFF_HemolysisSource = sim.getFieldFunctionManager().createFieldFunction();
        UserFieldFunction uFF_previousHemolysisSource = sim.getFieldFunctionManager().createFieldFunction();
        UserFieldFunction uFF_rotationSpeed = sim.getFieldFunctionManager().createFieldFunction();
        UserFieldFunction uFF_thresholdField = sim.getFieldFunctionManager().createFieldFunction();

        uFF_GWparam.getTypeOption().setSelected(FieldFunctionTypeOption.VECTOR);
        uFF_GWparam.setDefinition("[0.785,2.416,3.62e-7]");
        uFF_GWparam.setFunctionName("paramGW");
        uFF_GWparam.setPresentationName("Constant: GW Parameter");

        uFF_massFlow.getTypeOption().setSelected(FieldFunctionTypeOption.SCALAR);
        uFF_massFlow.setDefinition("0.01065");
        uFF_massFlow.setFunctionName("mflowRateInlet");
        uFF_massFlow.setPresentationName("Constant: Real MassFlowRate");

        uFF_thresholdValue.getTypeOption().setSelected(FieldFunctionTypeOption.SCALAR);
        uFF_thresholdValue.setDefinition("150");
        uFF_thresholdValue.setFunctionName("threshold");
        uFF_thresholdValue.setPresentationName("Constant: Threshold Stress");

        uFF_equivStress.getTypeOption().setSelected(FieldFunctionTypeOption.SCALAR);
        uFF_equivStress.setDefinition("${StrainRate}*${DynamicViscosity}");
        uFF_equivStress.setFunctionName("ssequiv");
        uFF_equivStress.setPresentationName("Equivalent Stress: Shear Scalar");

        uFF_HemolysisSource.getTypeOption().setSelected(FieldFunctionTypeOption.SCALAR);
        uFF_HemolysisSource.setDefinition("abs(${threshold} * $$paramGW[0] * pow($$paramGW[2],1/$$paramGW[0]) * pow($ssequiv, $$paramGW[1]/$$paramGW[0]) * pow(abs(${previousHemolysis} <= 1e-15?1e-15:${previousHemolysis})/${mflowRateInlet},1-1/$$paramGW[0]))");
        uFF_HemolysisSource.setFunctionName("hemolysisSource");
        uFF_HemolysisSource.setPresentationName("Hemolysis Source Term");

        uFF_previousHemolysisSource.getTypeOption().setSelected(FieldFunctionTypeOption.SCALAR);
        uFF_previousHemolysisSource.setDefinition("alternateValue(max(1e-15,abs(interpolatePositionTable(@Table(\"HemolysisPreviousStep\"),\"Hemolysis Source Term\"))), 1e-13)");
        uFF_previousHemolysisSource.setFunctionName("previousHemolysis");
        uFF_previousHemolysisSource.setPresentationName("previousHemolysis");

        uFF_rotationSpeed.getTypeOption().setSelected(FieldFunctionTypeOption.SCALAR);
        uFF_rotationSpeed.setDefinition("100");
        uFF_rotationSpeed.setFunctionName("rotationSpeed");
        uFF_rotationSpeed.setPresentationName("rotationSpeed");
        wallRelativeRotationProfile_disks.getMethod(FunctionScalarProfileMethod.class).setFieldFunction(uFF_rotationSpeed);

        uFF_thresholdField.getTypeOption().setSelected(FieldFunctionTypeOption.SCALAR);
        uFF_thresholdField.setDefinition("${ssequiv} > ${threshold} ?1:0");
        uFF_thresholdField.setFunctionName("thresholdss");
        uFF_thresholdField.setPresentationName("Threshold Field");

        // Create Table 
        Boundary boundary_inlet = region_fluid.getBoundaryManager().getBoundary("INLET");
        MassFlowRateProfile massFlowRateProfile_0 = boundary_inlet.getValues().get(MassFlowRateProfile.class);
        massFlowRateProfile_0.getMethod(ConstantScalarProfileMethod.class).getQuantity().setValue(0.01065);

        Boundary boundary_outlet = region_fluid.getBoundaryManager().getBoundary("OUTLET");
        // InterfaceBoundary interfaceBoundary_per0 = ((InterfaceBoundary) region_fluid.getBoundaryManager().getBoundary("SYM1 [Interface 1]"));
        // InterfaceBoundary interfaceBoundary_per1 = ((InterfaceBoundary) region_fluid.getBoundaryManager().getBoundary("SYM2 [Interface 1]"));
        Boundary boundary_volute = region_fluid.getBoundaryManager().getBoundary("VOLUTE");
        Boundary boundary_diff = region_fluid.getBoundaryManager().getBoundary("DIFFUSOR");

        XyzInternalTable xyzInternalTable_0 = sim.getTableManager().createTable(XyzInternalTable.class);
        xyzInternalTable_0.setPresentationName("HemolysisPreviousStep");
        xyzInternalTable_0.setFieldFunctions(new NeoObjectVector(new Object[]{uFF_HemolysisSource}));

        // xyzInternalTable_0.getParts().setObjects(region_fluid, boundary_per0, boundary_per1, boundary_inlet, boundary_outlet, interfaceBoundary_per0, interfaceBoundary_per1, boundary_volute, boundary_diff, boundary_disks);
        xyzInternalTable_0.getParts().setObjects(region_fluid, boundary_per0, boundary_per1, boundary_inlet, boundary_outlet, boundary_volute, boundary_diff, boundary_disks);

        TableUpdate tableUpdate_0 = xyzInternalTable_0.getTableUpdate();
        tableUpdate_0.getUpdateModeOption().setSelected(StarUpdateModeOption.ITERATION);
        IterationUpdateFrequency iterationUpdateFrequency_0 = tableUpdate_0.getIterationUpdateFrequency();
        iterationUpdateFrequency_0.setIterations(10);
        tableUpdate_0.setAutoExtract(true);

        UserFieldFunction uFF_freeEdgeStream = sim.getFieldFunctionManager().createFieldFunction();
        uFF_freeEdgeStream.getTypeOption().setSelected(FieldFunctionTypeOption.SCALAR);
        uFF_freeEdgeStream.setPresentationName("FreeStreamEdge");
        uFF_freeEdgeStream.setFunctionName("$WallDistance > 0.0001?1:0");

        GammaReThetaTransitionModel gammaReThetaTransitionModel_0 = physicsContinuum_fluid.getModelManager().getModel(GammaReThetaTransitionModel.class);
        gammaReThetaTransitionModel_0.setFreeStreamEdgeDefinitionFieldFunction(uFF_freeEdgeStream);

        // // Setting Source Term
        // region_fluid.getConditions().get(PassiveScalarUserSourceOption.class).setSelected(PassiveScalarUserSourceOption.INFERRED_DENSITY);
        // PassiveScalarUserSourceInferredDensity passiveScalarUserSourceInferredDensity_0 = region_fluid.getValues().get(PassiveScalarUserSourceInferredDensity.class);
        // passiveScalarUserSourceInferredDensity_0.setMethod(CompositeArrayProfileMethod.class);
        // ScalarProfile scalarProfile_blood = passiveScalarUserSourceInferredDensity_0.getMethod(CompositeArrayProfileMethod.class).getProfile(0);
        // scalarProfile_blood.setMethod(FunctionScalarProfileMethod.class);
        // scalarProfile_blood.getMethod(FunctionScalarProfileMethod.class).setFieldFunction(uFF_HemolysisSource);
        region_fluid.getConditions().get(PassiveScalarUserSourceOption.class).setSelected(PassiveScalarUserSourceOption.EXPLICIT_DENSITY);
        PassiveScalarUserSource passiveScalarUserSource_0 = region_fluid.getValues().get(PassiveScalarUserSource.class);
        passiveScalarUserSource_0.setMethod(CompositeArrayProfileMethod.class);

        ScalarProfile scalarProfile_0 = passiveScalarUserSource_0.getMethod(CompositeArrayProfileMethod.class).getProfile(0);
        scalarProfile_0.setMethod(FunctionScalarProfileMethod.class);
        scalarProfile_0.getMethod(FunctionScalarProfileMethod.class).setFieldFunction(uFF_HemolysisSource);

        // Make Reports, Monitors
        // PressureDrop
        PressureDropReport pressureDropReport = sim.getReportManager().createReport(PressureDropReport.class);
        pressureDropReport.getParts().setObjects(boundary_outlet);
        pressureDropReport.getLowPressureParts().setObjects(boundary_inlet);
        pressureDropReport.setPresentationName("presDrop");

        sim.getMonitorManager().createMonitorAndPlot(new NeoObjectVector(new Object[]{pressureDropReport}), true, "%1$s Plot");
        ReportMonitor reportMonitor_presDrop = ((ReportMonitor) sim.getMonitorManager().getMonitor("presDrop Monitor"));
        MonitorPlot monitorPlot_presDrop = sim.getPlotManager().createMonitorPlot(new NeoObjectVector(new Object[]{reportMonitor_presDrop}), "presDrop Monitor Plot");

        // Torque
        MomentReport momentReport_0 = sim.getReportManager().createReport(MomentReport.class);
        momentReport_0.setPresentationName("Torque");
        momentReport_0.getDirection().setComponents(0.0, 0.0, 1.0);
        momentReport_0.getParts().setObjects(boundary_disks);

        sim.getMonitorManager().createMonitorAndPlot(new NeoObjectVector(new Object[]{momentReport_0}), true, "%1$s Plot");
        ReportMonitor reportMonitor_Torque = ((ReportMonitor) sim.getMonitorManager().getMonitor("Torque Monitor"));
        MonitorPlot monitorPlot_Torque = sim.getPlotManager().createMonitorPlot(new NeoObjectVector(new Object[]{reportMonitor_Torque}), "Torque Monitor Plot");

        // Speed
        ExpressionReport expressionReport_rotationSpeed = sim.getReportManager().createReport(ExpressionReport.class);
        expressionReport_rotationSpeed.setDefinition("$rotationSpeed");
        expressionReport_rotationSpeed.setPresentationName("rotationSpeed");

        sim.getMonitorManager().createMonitorAndPlot(new NeoObjectVector(new Object[]{expressionReport_rotationSpeed}), true, "%1$s Plot");
        ReportMonitor reportMonitor_rotationSpeed = ((ReportMonitor) sim.getMonitorManager().getMonitor("rotationSpeed Monitor"));
        MonitorPlot monitorPlot_rotationSpeed = sim.getPlotManager().createMonitorPlot(new NeoObjectVector(new Object[]{reportMonitor_rotationSpeed}), "rotationSpeed Monitor Plot");

        // Hemolysis at Outlet
        PrimitiveFieldFunction primitiveFieldFunction_hemoRatio = ((PrimitiveFieldFunction) sim.getFieldFunctionManager().getFunction("HemolysisRatio"));

        MassFlowAverageReport massFlowAverageReport_hemolysis = sim.getReportManager().createReport(MassFlowAverageReport.class);
        massFlowAverageReport_hemolysis.setScalar(primitiveFieldFunction_hemoRatio);
        massFlowAverageReport_hemolysis.getParts().setObjects(boundary_outlet);
        massFlowAverageReport_hemolysis.setPresentationName("HIOutlet");

        sim.getMonitorManager().createMonitorAndPlot(new NeoObjectVector(new Object[]{massFlowAverageReport_hemolysis}), true, "%1$s Plot");
        ReportMonitor reportMonitor_hemolysis = ((ReportMonitor) sim.getMonitorManager().getMonitor("HIOutlet Monitor"));
        MonitorPlot monitorPlot_hemolysis = sim.getPlotManager().createMonitorPlot(new NeoObjectVector(new Object[]{reportMonitor_hemolysis}), "HIOutlet Monitor Plot");

        // Efficiency
        MaxReport maxReport_density = sim.getReportManager().createReport(MaxReport.class);
        maxReport_density.setPresentationName("density");

        PrimitiveFieldFunction ff_Density = ((PrimitiveFieldFunction) sim.getFieldFunctionManager().getFunction("Density"));
        maxReport_density.setScalar(ff_Density);
        maxReport_density.getParts().setObjects(boundary_inlet);

        ExpressionReport expressionReport_efficiency = sim.getReportManager().createReport(ExpressionReport.class);
        expressionReport_efficiency.setDefinition("-1*$presDropReport/$densityReport*$mflowRateInlet/($TorqueReport*$rotationSpeed)");
        expressionReport_efficiency.setPresentationName("Efficiency");

        sim.getMonitorManager().createMonitorAndPlot(new NeoObjectVector(new Object[]{expressionReport_efficiency}), true, "%1$s Plot");
        ReportMonitor reportMonitor_effciency = ((ReportMonitor) sim.getMonitorManager().getMonitor("Efficiency Monitor"));
        MonitorPlot monitorPlot_effciency = sim.getPlotManager().createMonitorPlot(new NeoObjectVector(new Object[]{reportMonitor_effciency}), "Efficiency Monitor Plot");

    }

    private void solve() {

        // get active session		
        Simulation sim = getActiveSimulation();

        // Some Constants
        double targetPressure = 13332.23;

        //some Variables
        double rotationRate = 2000;
        double maxJump = 50;
        double currentPressure;
        double jump;
        int nextIter = 10;

        boolean presOK = false;

        String tempFF;

        // Initiate Objects
        // Get Report pressure Drop			
        Report presDropReport = sim.getReportManager().getReport("presDrop");

        // Hemolysis Table for Previous step
        XyzInternalTable xyzInternalTable_0 = ((XyzInternalTable) sim.getTableManager().getTable("HemolysisPreviousStep"));
        TableUpdate tableUpdate_0 = xyzInternalTable_0.getTableUpdate();

        // // Get Wall Rotation Object
        // Region fluidregion = sim.getRegionManager().getRegion("FLUID");		
        // Boundary boundary_disks = fluidregion.getBoundaryManager().getBoundary("DISKS");		
        // WallRotationProfile wallRotationProfile = boundary_disks.getValues().get(WallRotationProfile.class);
        // get Solver Settings
        //get Solver settings
        StepStoppingCriterion maxStepStoppingCriterion
                = ((StepStoppingCriterion) sim.getSolverStoppingCriterionManager().getSolverStoppingCriterion("Maximum Steps"));
        MonitorIterationStoppingCriterion continuityStoppingCriterion
                = ((MonitorIterationStoppingCriterion) sim.getSolverStoppingCriterionManager().getSolverStoppingCriterion("Continuity Criterion"));
        MonitorIterationStoppingCriterion xMomentumStoppingCriterion
                = ((MonitorIterationStoppingCriterion) sim.getSolverStoppingCriterionManager().getSolverStoppingCriterion("X-momentum Criterion"));
        MonitorIterationStoppingCriterion yMomentumStoppingCriterion
                = ((MonitorIterationStoppingCriterion) sim.getSolverStoppingCriterionManager().getSolverStoppingCriterion("Y-momentum Criterion"));
        MonitorIterationStoppingCriterion zMomentumStoppingCriterion
                = ((MonitorIterationStoppingCriterion) sim.getSolverStoppingCriterionManager().getSolverStoppingCriterion("Z-momentum Criterion"));
        MonitorIterationStoppingCriterion tkeStoppingCriterion
                = ((MonitorIterationStoppingCriterion) sim.getSolverStoppingCriterionManager().getSolverStoppingCriterion("Tke Criterion"));
        MonitorIterationStoppingCriterion sdrStoppingCriterion
                = ((MonitorIterationStoppingCriterion) sim.getSolverStoppingCriterionManager().getSolverStoppingCriterion("Sdr Criterion"));
        MonitorIterationStoppingCriterion rthetaStoppingCriterion
                = ((MonitorIterationStoppingCriterion) sim.getSolverStoppingCriterionManager().getSolverStoppingCriterion("ReTheta_t Criterion"));
        MonitorIterationStoppingCriterion intermittStoppingCriterion
                = ((MonitorIterationStoppingCriterion) sim.getSolverStoppingCriterionManager().getSolverStoppingCriterion("Intermittency Criterion"));

        UserFieldFunction currentSpeedFF = ((UserFieldFunction) sim.getFieldFunctionManager().getFunction("rotationSpeed"));

        UserFieldFunction hemolysisPrevStepFF = ((UserFieldFunction) sim.getFieldFunctionManager().getFunction("previousHemolysis"));
        tempFF = hemolysisPrevStepFF.getDefinition();

        // Passive Scalar Solver
        PassiveScalarSolver passiveScalarSolver_0 = ((PassiveScalarSolver) sim.getSolverManager().getSolver(PassiveScalarSolver.class));
        passiveScalarSolver_0.setFrozen(true);

        WallDistanceSolver wallDistanceSolver_0 = ((WallDistanceSolver) sim.getSolverManager().getSolver(WallDistanceSolver.class));

        // Init
        tableUpdate_0.setAutoExtract(false);
        hemolysisPrevStepFF.setDefinition("$thresholdss*1e-15");
        xyzInternalTable_0.extract();
        sim.getSimulationIterator().step(10);
        passiveScalarSolver_0.setFrozen(false);
        sim.getSimulationIterator().step(10);
        xyzInternalTable_0.extract();
        hemolysisPrevStepFF.setDefinition(tempFF);
        sim.getSimulationIterator().step(10);
        xyzInternalTable_0.extract();
        tableUpdate_0.setAutoExtract(true);
        sim.getSimulationIterator().step(10);
        wallDistanceSolver_0.setFrozen(true);
        sim.getSimulationIterator().step(nextIter);

        // Solve
        do {

            // Get Pressuredrop
            currentPressure = presDropReport.getReportMonitorValue();

            // Get current Rotation speed
            rotationRate = Float.valueOf(currentSpeedFF.getDefinition());
            //rotationRate = wallRotationProfile.getMethod(ConstantScalarProfileMethod.class).getQuantity().getValue();

            // Calc new Roation Speed
            if (currentPressure > 0) {
                jump = Math.sqrt(targetPressure / currentPressure) * rotationRate - rotationRate;
                if (Math.pow(jump, 2) > Math.pow(maxJump, 2)) {
                    rotationRate = rotationRate + Math.copySign(maxJump, jump);
                } else {
                    rotationRate = rotationRate + jump * 0.05;
                }
                nextIter = 100;
            } else {
                // If Turbine effect
                rotationRate = rotationRate + 10.0;
                nextIter = 100;

            }

            // set current speed
            //wallRotationProfile.getMethod(ConstantScalarProfileMethod.class).getQuantity().setValue(rotationRate);
            currentSpeedFF.setDefinition(String.valueOf(rotationRate));

            // Advance solution
            sim.getSimulationIterator().step(nextIter);

            // Get Pressuredrop
            currentPressure = presDropReport.getReportMonitorValue();

            if (currentPressure < (targetPressure + 10.0) && currentPressure > (targetPressure - 10.0)) {
                presOK = true;
            } else {
                presOK = false;
            }
	// do as long as 
        } while (!maxStepStoppingCriterion.getIsSatisfied() && !(presOK && continuityStoppingCriterion.getIsSatisfied() && xMomentumStoppingCriterion.getIsSatisfied() && yMomentumStoppingCriterion.getIsSatisfied() && zMomentumStoppingCriterion.getIsSatisfied() && tkeStoppingCriterion.getIsSatisfied() && sdrStoppingCriterion.getIsSatisfied()));

        // Disable Stopping Criteria
        StepStoppingCriterion stepStoppingCriterion_0 = ((StepStoppingCriterion) sim.getSolverStoppingCriterionManager().getSolverStoppingCriterion("Maximum Steps"));
        stepStoppingCriterion_0.setIsUsed(false);

        sim.getSimulationIterator().step(2000);

    }

    private void save() {

        Simulation sim = getActiveSimulation();

        // Save Simulation data
        sim.saveState(resolvePath("slice.sim"));

        // Save Screens	
        MonitorPlot monitorPlot_0 = ((MonitorPlot) sim.getPlotManager().getPlot("presDrop Monitor Plot"));
        MonitorPlot monitorPlot_1 = ((MonitorPlot) sim.getPlotManager().getPlot("Torque Monitor Plot"));
        MonitorPlot monitorPlot_2 = ((MonitorPlot) sim.getPlotManager().getPlot("rotationSpeed Monitor Plot"));
        MonitorPlot monitorPlot_3 = ((MonitorPlot) sim.getPlotManager().getPlot("HIOutlet Monitor Plot"));
        MonitorPlot monitorPlot_4 = ((MonitorPlot) sim.getPlotManager().getPlot("Efficiency Monitor Plot"));
        ResidualPlot residualPlot_0 = ((ResidualPlot) sim.getPlotManager().getPlot("Residuals"));

        monitorPlot_0.export(resolvePath("presDrop.csv"), ",");
        monitorPlot_1.export(resolvePath("torque.csv"), ",");
        monitorPlot_2.export(resolvePath("rotationSpeed.csv"), ",");
        monitorPlot_3.export(resolvePath("HIoutlet.csv"), ",");
        monitorPlot_4.export(resolvePath("efficiency.csv"), ",");

        monitorPlot_0.encode(resolvePath("presDrop.png"), "png", 1280, 1024);
        monitorPlot_1.encode(resolvePath("torque.png"), "png", 1280, 1024);
        monitorPlot_2.encode(resolvePath("rotationSpeed.png"), "png", 1280, 1024);
        monitorPlot_3.encode(resolvePath("HIoutlet.png"), "png", 1280, 1024);
        monitorPlot_4.encode(resolvePath("efficiency.png"), "png", 1280, 1024);
        residualPlot_0.encode(resolvePath("residuals.png"), "png", 1280, 1024);

    }

}
