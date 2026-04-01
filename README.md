# Partial Equilibrium Model of New and Used Cars

A dynamic partial-equilibrium model of the car market implemented in Julia using [SquareModels](https://github.com/MartinBonde/SquareModels) and [Ipopt](https://github.com/coin-or/Ipopt). The model captures substitution between new and used cars across fuel types (petrol and electric), with habit formation in the used-car market.

## Model Structure

The household allocates total consumption $`C_t`$ between a car-service aggregate $`d_t`$ and non-car consumption $`c_{nc,t}`$ via a nested CES demand system:

```
Total consumption C
├── Non-car consumption (c_nc)
└── Car-service aggregate (d)
    ├── New cars (d_new)
    │   ├── Petrol
    │   └── Electric
    └── Used cars (d_used)
        ├── Petrol (with habit)
        └── Electric (with habit)
```

Total consumption $`C_t`$ and the non-car price $`p_{nc,t}`$ are exogenous. The model determines the car/non-car split, allocation across new/used and fuel types, and the market-clearing prices.

## Utility Function

The representative household maximises

$$\sum_{t=0}^{\infty} \beta^t u(C_t)$$

where period utility is a nested CES aggregate of total consumption $`C_t`$:

$$C_t = \left[\mu_{nc}^{1/\sigma_C} c_{nc,t}^{(\sigma_C-1)/\sigma_C} + \mu_d^{1/\sigma_C} d_t^{(\sigma_C-1)/\sigma_C}\right]^{\sigma_C/(\sigma_C-1)}$$

The car-service aggregate $`d_t`$ nests new and used cars:

$$d_t = \left[(\mu^{new})^{1/\sigma} (d_t^{new})^{(\sigma-1)/\sigma} + (\mu^{used})^{1/\sigma} (d_t^{used})^{(\sigma-1)/\sigma}\right]^{\sigma/(\sigma-1)}$$

New cars aggregate over fuel types:

$$d_t^{new} = \left[\sum_f (\mu_f^{new})^{1/\sigma^{new}} (d_{f,t}^{new})^{(\sigma^{new}-1)/\sigma^{new}}\right]^{\sigma^{new}/(\sigma^{new}-1)}$$

Used cars aggregate over fuel types with **habit formation**:

$$d_t^{used} = \left[\sum_f (\mu_f^{used})^{1/\sigma^{used}} \left(d_{f,t}^{used} - h_f d_{f,t-1}^{used}\right)^{(\sigma^{used}-1)/\sigma^{used}}\right]^{\sigma^{used}/(\sigma^{used}-1)}$$

The habit term $`h_f d_{f,t-1}^{used}`$ is the reference point: only the stock of fuel type $`f`$ in excess of this generates marginal utility. This creates inertia in the fuel-type composition of the used fleet — a household inheriting a large petrol stock finds it costly to shrink it because the habit-adjusted quantity $`d_{f,t}^{used} - h_f d_{f,t-1}^{used}`$ falls, depressing utility.

## Stock Accumulation

Let $`a = 0`$ denote new cars so that $`d_{0,f,t} \equiv d^{new}_{f,t}`$. Cars of age $`a`$ depreciate at rate $`\delta_{a,f,t}`$ (which can vary by age, fuel type, and time) and the stock evolves as:

$$d_{a,f,t} = (1 - \delta_{a-1,f,t}) d_{a-1,f,t-1}, \quad a \geq 1$$

Used cars of a given fuel type are perfect substitutes, $`d^{used}_{f,t} = \sum_{a=1}^{\infty} d_{a,f,t}`$. If depreciation rates are age-independent for $`a \geq 1`$ (i.e. $`\delta_{a,f,t} = \delta_{f,t}`$ for all $`a \geq 1`$, while $`\delta_{0,f,t}`$ may differ), the stock simplifies to:

$$d^{used}_{f,t} = (1 - \delta_{f,t}) d^{used}_{f,t-1} + (1 - \delta_{0,f,t}) d^{new}_{f,t-1}$$

## Demand System

### Top level: cars vs. non-car consumption

$$d_t = \mu_d C_t \left(\frac{p^d_t}{p^C_t}\right)^{-\sigma_C}, \qquad c_{nc,t} = \mu_{nc} C_t \left(\frac{p_{nc,t}}{p^C_t}\right)^{-\sigma_C}$$

$$p^C_t C_t = p^d_t d_t + p_{nc,t} c_{nc,t}$$

### Car nest: new vs. used

$$d^{new}_t = \mu^{new} d_t \left(\frac{p^{uc,new}_t}{p^d_t}\right)^{-\sigma}, \qquad d^{used}_t = \mu^{used} d_t \left(\frac{p^{uc,used}_t}{p^d_t}\right)^{-\sigma}$$

$$p^d_t d_t = p^{uc,new}_t d^{new}_t + p^{uc,used}_t d^{used}_t$$

### New-car nest: across fuel types

$$d^{new}_{f,t} = \mu^{new}_f d^{new}_t \left(\frac{p^{uc,new}_{f,t}}{p^{uc,new}_t}\right)^{-\sigma^{new}}$$

$$p^{uc,new}_t d^{new}_t = \sum_f p^{uc,new}_{f,t} d^{new}_{f,t}$$

The user cost of a new car of fuel type $`f`$ is:

$$p^{uc,new}_{f,t} = p^{new}_{f,t} - \frac{1 - \delta_{0,f,t+1}}{1 + r_{t+1}} p^{used}_{f,t+1}$$

### Used-car nest: across fuel types (with habits)

The CES aggregator in this nest is defined over habit-adjusted quantities $`d^{used}_{f,t} - h_f d^{used}_{f,t-1}`$, where $`h_f \in (0,1)`$ is a habit parameter. Only the stock in excess of the habit reference point generates marginal utility, creating inertia: a household with a large inherited stock of fuel type $`f`$ finds it costly to reduce holdings.

$$d^{used}_{f,t} - h_f d^{used}_{f,t-1} = \mu^{used}_f d^{used}_t \left(\frac{p^{uc}_{f,t}}{p^{uc,used}_t}\right)^{-\sigma^{used}}$$

$$p^{uc,used}_t d^{used}_t = \sum_f p^{uc}_{f,t} \left(d^{used}_{f,t} - h_f d^{used}_{f,t-1}\right)$$

The user cost of a used car of fuel type $`f`$ is:

$$p^{uc}_{f,t} = p^{used}_{f,t} - \frac{1 - \delta_{f,t+1}}{1 + r_{t+1}} p^{used}_{f,t+1} + \frac{1 - \delta_{f,t+1}}{1 + r_{t+1}} h_f p^{uc,used}_{t+1} \mu^{used}_f \left(\frac{d^{used}_{t+1}}{d^{used}_{f,t+1} - h_f d^{used}_{f,t}}\right)^{1/\sigma^{used}}$$

The third term is the **habit premium**: holding more used cars of type $`f`$ today raises next period's reference point by $`h_f(1 - \delta_{f,t+1})`$, reducing the effective service flow and increasing the marginal cost of maintaining the same utility level tomorrow. When $`h_f = 0`$ the habit premium vanishes.

## Calibration

The share parameters $`\mu_d, \mu_{nc}, \mu^{new}, \mu^{used}, \mu^{new}_f, \mu^{used}_f`$ are calibrated by swapping them for initial-period quantities and solving for the values that match base-year data.

### Baseline Parameters

| Parameter | Value | Description |
|-----------|-------|-------------|
| $`\sigma_C`$ | 0.5 | Elasticity: cars vs. non-car |
| $`\sigma`$ | 1.5 | Elasticity: new vs. used |
| $`\sigma^{new}`$ | 3.0 | Elasticity: across fuel types (new) |
| $`\sigma^{used}`$ | 3.0 | Elasticity: across fuel types (used) |
| $`h_f`$ | 0.8 | Habit parameter (both fuel types) |
| $`r`$ | 0.04 | Interest rate |
| $`\delta_0`$ | 0.25 | First-period depreciation (new to used) |
| $`\delta`$ | 0.10 | Ongoing used-car depreciation |

## Scenario 1: EV Subsidy

The first scenario simulates a **10% reduction in the purchase price of electric cars** from 2026 onward. The unfinanced subsidy shifts new-car demand toward electric vehicles, gradually building up the electric used-car stock as cheaper EVs flow through the depreciation pipeline, while all equilibrium prices adjust.

### Results

![EV Subsidy Scenario](cars_scenario1.svg)

- **New cars by fuel type** — Electric purchases rise while petrol purchases contract.
- **Used-car stock by fuel type** — The stock responds with a lag as new electric cars depreciate into the used market; the petrol stock declines as fewer petrol cars enter the pipeline.
- **Aggregates** — Total new-car purchases increase, funded by the unfinanced subsidy.
- **User costs by fuel type** — The user cost of new electric cars drops directly with the price reduction; used-car user costs adjust endogenously as stock composition shifts.
- **Used-car spot prices by fuel type** — The growing stock of electric vehicles depresses electric resale values; used petrol prices edge up as the petrol stock thins.
- **New-car purchase prices by fuel type** — The exogenous shock: a flat 10% reduction for electric cars, petrol unchanged.

## Scenario 2: PV-Neutral Tax/Subsidy

The second scenario pairs the **10% EV subsidy** with a **constant endogenous ad-valorem tax on petrol cars** $`\tau_{petrol}`$, calibrated so that the present value of net tax revenue is exactly zero:

$$\sum_t \frac{1}{(1+r_t)^{t-t_1}} \sum_f \tau_{f,t} p^{new}_{f,t} d^{new}_{f,t} = 0$$

The consumer-facing purchase price becomes $`p^{new}_{f,t}(1 + \tau_{f,t})`$, which enters the user cost of new cars.

The required petrol tax rate depends critically on the substitution elasticities. High fuel-type substitutability ($`\sigma^{new} = 3`$) means the EV subsidy erodes the petrol tax base aggressively — households switch away from petrol cars easily, shrinking the revenue that any given tax rate can raise. The petrol tax must therefore be *higher* than the 10% subsidy it finances. More generally, $`\tau_{petrol}`$ is an increasing function of $`\sigma^{new}`$: the easier it is to substitute between fuel types, the more the tax base shrinks, and the higher the rate needed to close the budget. In the limit $`\sigma^{new} \to \infty`$, the tax base vanishes entirely and no finite rate can balance the budget.

### Results

![PV-Neutral Scenario](cars_scenario2.svg)

- **New cars by fuel type** — The simultaneous tax on petrol and subsidy on electric drives a larger compositional shift than Scenario 1, since both margins push in the same direction.
- **Used-car stock by fuel type** — The petrol used-car stock declines more steeply as the tax chokes off new petrol inflows at the source.
- **Aggregates** — Despite fiscal neutrality, the car aggregate $`d_t`$ falls. Both the subsidy and the tax distort relative prices away from their undistorted values, and the CES aggregator registers the combined efficiency loss as a decline in the car-service aggregate. Revenue neutrality balances the budget, not welfare.
- **User costs by fuel type** — A symmetric wedge opens up: petrol user costs rise due to the tax while electric user costs fall, creating a wider gap than the subsidy-only scenario.
- **Used-car spot prices by fuel type** — The petrol resale price *rises* as reduced future supply makes the surviving stock more scarce, while the electric resale price falls as subsidized vehicles flood the secondary market.
- **Implied tax/subsidy rates** — The constant −10% electric subsidy and the endogenous petrol tax that balances revenue in present value.

## Running

```julia
julia --project=. cars.jl
```

This solves the baseline calibration, runs both counterfactual scenarios, and saves plots to `cars_scenario1.svg` and `cars_scenario2.svg`.

### Dependencies

- [JuMP](https://github.com/jump-dev/JuMP.jl) + [Ipopt](https://github.com/jump-dev/Ipopt.jl) — nonlinear optimization
- [SquareModels](https://github.com/MartinBonde/SquareModels) — model definition and solution framework
- [CairoMakie](https://github.com/MakieOrg/Makie.jl) — plotting
