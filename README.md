# Partial Equilibrium Model of New and Used Cars

A dynamic partial-equilibrium model of the car market implemented in Julia using [SquareModels](https://github.com/MartinBonde/SquareModels) and [Ipopt](https://github.com/coin-or/Ipopt). The model captures substitution between new and used cars across fuel types (petrol and electric), with habit formation in the used-car market.

## Model Structure

The household allocates total consumption \( C_t \) between a car-service aggregate \( d_t \) and non-car consumption \( c_{nc,t} \) via a CES demand system. The car aggregate is further decomposed through a nested CES structure:

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

Total consumption \( C_t \) and the non-car price \( p_{nc,t} \) are exogenous. The model determines the car/non-car split, allocation across new/used and fuel types, and the market-clearing prices.

## Demand System

### Top level: cars vs. non-car consumption

The CES demands for the car aggregate and non-car consumption are:

$$d_t = \mu_d \, C_t \left(\frac{p^{uc,d}_t}{p^C_t}\right)^{-\sigma_C}, \qquad c_{nc,t} = \mu_{nc} \, C_t \left(\frac{p_{nc,t}}{p^C_t}\right)^{-\sigma_C}$$

subject to the budget constraint:

$$p^C_t \, C_t = p^{uc,d}_t \, d_t + p_{nc,t} \, c_{nc,t}$$

### Car nest: new vs. used

Let \( \sigma \) be the elasticity of substitution between new and used cars, with share parameters \( \mu^{new} \) and \( \mu^{used} \). The conditional demands are:

$$d^{new}_t = \mu^{new} \, d_t \left(\frac{p^{uc,new}_t}{p^d_t}\right)^{-\sigma}, \qquad d^{used}_t = \mu^{used} \, d_t \left(\frac{p^{uc,used}_t}{p^d_t}\right)^{-\sigma}$$

with the budget constraint:

$$p^d_t \, d_t = p^{uc,new}_t \, d^{new}_t + p^{uc,used}_t \, d^{used}_t$$

The user cost of the car-service aggregate is:

$$p^{uc,d}_t = p^d_t - \frac{1 - \delta^d_{t+1}}{1 + r_{t+1}} \, p^d_{t+1}$$

### New-car nest: across fuel types

Let \( \sigma^{new} \) be the elasticity and \( \mu^{new}_f \) the share parameters. The conditional demands are:

$$d^{new}_{f,t} = \mu^{new}_f \, d^{new}_t \left(\frac{p^{uc,new}_{f,t}}{p^{uc,new}_t}\right)^{-\sigma^{new}}$$

The user cost of a new car of fuel type \( f \) is:

$$p^{uc,new}_{f,t} = p^{new}_{f,t} - \frac{1 - \delta_{0,f,t+1}}{1 + r_{t+1}} \, p^{used}_{f,t+1}$$

A fraction \( 1 - \delta_{0,f,t+1} \) of each new car survives into the used stock next period — the same first-period depreciation rate that enters the accumulation equation — so the resale value is discounted accordingly.

### Used-car nest: across fuel types (with habits)

The used-car nest introduces habit formation. The CES aggregator is defined over habit-adjusted quantities \( d^{used}_{f,t} - h_f \, d^{used}_{f,t-1} \), where \( h_f \in (0,1) \) is a habit parameter. Only the stock in excess of the habit reference point generates marginal utility.

The conditional demands are:

$$d^{used}_{f,t} - h_f \, d^{used}_{f,t-1} = \mu^{used}_f \, d^{used}_t \left(\frac{p^{uc}_{f,t}}{p^{uc,used}_t}\right)^{-\sigma^{used}}$$

The budget constraint is:

$$p^{uc,used}_t \, d^{used}_t = \sum_f p^{uc}_{f,t} \left(d^{used}_{f,t} - h_f \, d^{used}_{f,t-1}\right)$$

The user cost of a used car of fuel type \( f \) is:

$$p^{uc}_{f,t} = p^{used}_{f,t} - \frac{1 - \delta_{f,t+1}}{1 + r_{t+1}} \, p^{used}_{f,t+1} + \frac{1 - \delta_{f,t+1}}{1 + r_{t+1}} \, h_f \, p^{uc,used}_{t+1} \, \mu^{used}_f \left(\frac{d^{used}_{t+1}}{d^{used}_{f,t+1} - h_f \, d^{used}_{f,t}}\right)^{1/\sigma^{used}}$$

The first two terms are the standard asset-pricing user cost: the purchase price minus the discounted resale value. The third term is the **habit premium**: holding more used cars of type \( f \) today raises next period's reference point by \( h_f(1 - \delta_{f,t+1}) \), which reduces the effective service flow \( d^{used}_{f,t+1} - h_f \, d^{used}_{f,t} \) and hence increases the marginal cost of maintaining the same utility level tomorrow. When \( h_f = 0 \) the habit premium vanishes and the user cost reduces to the standard expression.

## Stock Accumulation

Cars of age \( a \) depreciate at rate \( \delta_{a,f,t} \), which can vary by age, fuel type, and time. The stock evolves as:

$$d_{a,f,t} = (1 - \delta_{a-1,f,t}) \, d_{a-1,f,t-1}, \quad a \geq 1$$

Used cars of a given fuel type are perfect substitutes so that \( d^{used}_{f,t} = \sum_{a=1}^{\infty} d_{a,f,t} \). If depreciation rates are independent of age for \( a \geq 1 \) (i.e. \( \delta_{a,f,t} = \delta_{f,t} \) for \( a \geq 1 \), while \( \delta_{0,f,t} \) may differ), the stock of used cars simplifies to:

$$d^{used}_{f,t} = (1 - \delta_{f,t}) \, d^{used}_{f,t-1} + (1 - \delta_{0,f,t}) \, d^{new}_{f,t-1}$$

## Calibration

The model is calibrated by swapping share parameters for initial-period quantities. The share parameters \( \mu_d, \mu_{nc}, \mu^{new}, \mu^{used}, \mu^{new}_f, \mu^{used}_f \) are solved to match base-year data on car demands and prices.

### Baseline Parameters

| Parameter | Value | Description |
|-----------|-------|-------------|
| \( \sigma_C \) | 0.5 | Elasticity: cars vs. non-car |
| \( \sigma \) | 1.5 | Elasticity: new vs. used |
| \( \sigma^{new} \) | 3.0 | Elasticity: across fuel types (new) |
| \( \sigma^{used} \) | 3.0 | Elasticity: across fuel types (used) |
| \( h_f \) | 0.3 | Habit parameter (both fuel types) |
| \( r \) | 0.03 | Interest rate |
| \( \delta^d \) | 0.10 | Aggregate car-service depreciation |
| \( \delta_0 \) | 0.25 | First-period depreciation (new to used) |
| \( \delta \) | 0.10 | Ongoing used-car depreciation |

## Counterfactual: EV Subsidy

The scenario simulates a **20% reduction in the purchase price of electric cars** from 2026 to 2040. The model traces how this subsidy propagates through the nested demand system — shifting new-car purchases toward electric vehicles, gradually building up the electric used-car stock, and adjusting all equilibrium prices.

### Results

![EV Subsidy Scenario](cars_scenario.svg)

The figure shows four panels:

- **New cars by fuel type** — The electric car subsidy increases new electric car purchases while crowding out petrol cars.
- **Used-car stock by fuel type** — The expanded flow of new electric cars gradually raises the electric used-car stock; the petrol used-car stock declines as fewer petrol cars enter the pipeline.
- **Aggregates** — Total new-car demand rises modestly; the used-car stock shifts in composition. Non-car consumption adjusts via the top-level CES.
- **User costs by fuel type** — The user cost of new electric cars drops directly with the subsidy. Used-car prices adjust endogenously as the composition of the stock changes.

## Running

```julia
julia --project=. cars.jl
```

This solves the baseline calibration, runs the EV subsidy counterfactual, and saves the plot to `cars_scenario.svg`.

### Dependencies

- [JuMP](https://github.com/jump-dev/JuMP.jl) + [Ipopt](https://github.com/jump-dev/Ipopt.jl) — nonlinear optimization
- [SquareModels](https://github.com/MartinBonde/SquareModels) — model definition and solution framework
- [CairoMakie](https://github.com/MakieOrg/Makie.jl) — plotting
