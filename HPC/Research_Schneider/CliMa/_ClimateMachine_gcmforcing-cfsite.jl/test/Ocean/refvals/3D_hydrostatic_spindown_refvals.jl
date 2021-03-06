# [
#  [ MPIStateArray Name, Field Name, Maximum, Minimum, Mean, Standard Deviation ],
#  [         :                :          :        :      :          :           ],
# ]
parr = [
    ["Q", "u[1]", 12, 12, 0, 12],
    ["Q", "u[2]", 0, 0, 0, 0],
    ["Q", :η, 12, 12, 0, 12],
    ["Q", :θ, 15, 15, 15, 15],
    ["s_aux", :y, 15, 15, 15, 15],
    ["s_aux", :w, 12, 12, 0, 12],
    ["s_aux", :pkin, 15, 15, 15, 15],
    ["s_aux", :wz0, 12, 12, 0, 12],
    ["aux", "uᵈ[1]", 15, 15, 15, 15],
    ["aux", "uᵈ[2]", 15, 15, 15, 15],
    ["aux", "ΔGᵘ[1]", 15, 15, 15, 15],
    ["aux", "ΔGᵘ[2]", 15, 15, 15, 15],
]

### fully explicit
explicit = [
    [
        "state",
        "u[1]",
        -9.58544066049463849843e-01,
        9.58544066049465071089e-01,
        -6.13908923696726568442e-17,
        4.45400263687296238402e-01,
    ],
    [
        "state",
        "u[2]",
        -4.33260855197780914020e-14,
        1.99397359064900580713e-14,
        -1.30343101521973575132e-16,
        2.95525689323845820606e-15,
    ],
    [
        "state",
        :η,
        -8.52732886154656810618e-01,
        8.52845586939211197652e-01,
        2.20052243093959998331e-14,
        6.02992088522925295813e-01,
    ],
    [
        "state",
        :θ,
        0.00000000000000000000e+00,
        0.00000000000000000000e+00,
        0.00000000000000000000e+00,
        0.00000000000000000000e+00,
    ],
    [
        "aux",
        :y,
        0.00000000000000000000e+00,
        1.00000000000000011642e+06,
        5.00000000000000000000e+05,
        2.92775877460665535182e+05,
    ],
    [
        "aux",
        :w,
        -4.04553460063758398447e-04,
        4.04714358463272711169e-04,
        4.75730566051879549438e-19,
        1.63958655681888576441e-04,
    ],
    [
        "aux",
        :pkin,
        0.00000000000000000000e+00,
        0.00000000000000000000e+00,
        0.00000000000000000000e+00,
        0.00000000000000000000e+00,
    ],
    [
        "aux",
        :wz0,
        -2.01164684799271339293e-04,
        2.01041968159484089494e-04,
        -2.10942374678779754294e-20,
        1.42228420244455277133e-04,
    ],
    [
        "aux",
        "uᵈ[1]",
        0.00000000000000000000e+00,
        0.00000000000000000000e+00,
        0.00000000000000000000e+00,
        0.00000000000000000000e+00,
    ],
    [
        "aux",
        "uᵈ[2]",
        0.00000000000000000000e+00,
        0.00000000000000000000e+00,
        0.00000000000000000000e+00,
        0.00000000000000000000e+00,
    ],
    [
        "aux",
        "ΔGᵘ[1]",
        0.00000000000000000000e+00,
        0.00000000000000000000e+00,
        0.00000000000000000000e+00,
        0.00000000000000000000e+00,
    ],
    [
        "aux",
        "ΔGᵘ[2]",
        0.00000000000000000000e+00,
        0.00000000000000000000e+00,
        0.00000000000000000000e+00,
        0.00000000000000000000e+00,
    ],
]

imex = [
    [
        "state",
        "u[1]",
        -9.57705359685807500192e-01,
        9.57705359685806500991e-01,
        4.54747350886464140750e-17,
        4.45323727947873004851e-01,
    ],
    [
        "state",
        "u[2]",
        -4.07596731389788922865e-14,
        2.53992382675900717845e-14,
        9.88087105439151860723e-18,
        2.63742927063789745855e-15,
    ],
    [
        "state",
        :η,
        -8.55941716778505390373e-01,
        8.56014908798837681481e-01,
        2.16732587432488803528e-14,
        6.05260814296588400829e-01,
    ],
    [
        "state",
        :θ,
        0.00000000000000000000e+00,
        0.00000000000000000000e+00,
        0.00000000000000000000e+00,
        0.00000000000000000000e+00,
    ],
    [
        "aux",
        :y,
        0.00000000000000000000e+00,
        1.00000000000000011642e+06,
        5.00000000000000000000e+05,
        2.92775877460665535182e+05,
    ],
    [
        "aux",
        :w,
        -4.03419301824637895407e-04,
        4.03626959596159180614e-04,
        -3.33066907387546970565e-21,
        1.63816237485524198066e-04,
    ],
    [
        "aux",
        :pkin,
        0.00000000000000000000e+00,
        0.00000000000000000000e+00,
        0.00000000000000000000e+00,
        0.00000000000000000000e+00,
    ],
    [
        "aux",
        :wz0,
        -1.97016630922629948884e-04,
        1.97066086934914155761e-04,
        -7.49400541621980656312e-19,
        1.39393318574472925390e-04,
    ],
    [
        "aux",
        "uᵈ[1]",
        0.00000000000000000000e+00,
        0.00000000000000000000e+00,
        0.00000000000000000000e+00,
        0.00000000000000000000e+00,
    ],
    [
        "aux",
        "uᵈ[2]",
        0.00000000000000000000e+00,
        0.00000000000000000000e+00,
        0.00000000000000000000e+00,
        0.00000000000000000000e+00,
    ],
    [
        "aux",
        "ΔGᵘ[1]",
        0.00000000000000000000e+00,
        0.00000000000000000000e+00,
        0.00000000000000000000e+00,
        0.00000000000000000000e+00,
    ],
    [
        "aux",
        "ΔGᵘ[2]",
        0.00000000000000000000e+00,
        0.00000000000000000000e+00,
        0.00000000000000000000e+00,
        0.00000000000000000000e+00,
    ],
]

refVals = (explicit = (explicit, parr), imex = (imex, parr))
