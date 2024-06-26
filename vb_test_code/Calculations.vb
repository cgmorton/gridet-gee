Public Module Program

  Public Const π As Double = Math.PI

  'Derived from Strong, C., Khatri, K. B., Kochanski, A. K., Lewis, C. S., & Allen, L. N. (2017). Reference evapotranspiration from coarse-scale and dynamically downscaled data in complex terrain: Sensitivity to interpolation and resolution. Journal of Hydrology, 548, 406-418.
  Public MonthlyLapseRate() As Double = {-0.001742, -0.002699, -0.003064, -0.003432, -0.003262, -0.00319, -0.003046, -0.002941, -0.002659, -0.002622, -0.002247, -0.002132} 'Fahrenheit/Foot
  Public AirTemperatureCorrectionCoefficients As Double() = {1.58, 0.59, -1.53, -3.73, 1.4, 0.0551}
  Public HumidityCorrectionCoefficients As Double() = {-21.9, 0.78, 3.55, 11.6, -5.05, 0.274}
  Public NLDAS_2AStartDate As DateTime = New DateTime(1979, 1, 1, 13, 0, 0, 0, DateTimeKind.Utc)


  Function CalculateReferenceEvapotranspirationNLDAS_2A(
    ByVal MinDate As DateTime,
    ByVal MaxDate As DateTime,
    ByVal Longitude As Double,
    ByVal Latitude As Double,
    ByVal Elevation As Double,
    ByVal Slope As Double,
    ByVal Aspect As Double,
    ByVal AirTemperature As Double,
    ByVal AirPressure As Double,
    ByVal NLDAS_Longitude As Double,
    ByVal NLDAS_Latitude As Double,
    ByVal NLDAS_Elevation As Double,
    ByVal NLDAS_Slope As Double,
    ByVal NLDAS_Aspect As Double,
    ByVal NLDAS_AirTemperature As Double,
    ByVal NLDAS_AirPressure As Double,
    ByVal NLDAS_Rs As Double,
    ByVal NLDAS_RH As Double,
    ByVal NLDAS_WindU As Double,
    ByVal NLDAS_WindV As Double
  ) As Double

    'Determine Calculation Hours
    Dim First As Integer = MinDate.Subtract(NLDAS_2AStartDate).TotalHours
    Dim Last As Integer = First + (MaxDate.Subtract(MinDate).TotalHours \ 24) * 24

    Dim Hour As Integer = First

    Dim RecordDate = NLDAS_2AStartDate.AddHours(Hour).ToUniversalTime

    'Calculate Solar Positioning for Multiple Time Steps
    Dim TimeStep = 15 'minutes
    Dim HourlySteps As Integer = 60 / TimeStep
    Dim SolarPositions(1440 / TimeStep) As SolarPosition
    For I As Integer = 0 To SolarPositions.Length - 1
      SolarPositions(I) = CalculateSolarPosition(RecordDate.AddMinutes(I * TimeStep - 30))
    Next

    'Calculate NLDAS-2A Extraterrestrial Radiation
    Dim NLDAS_RaPixels(SolarPositions.Length - 1) As Single
    For J As Integer = 0 To SolarPositions.Length - 1
      Dim H_temp As Double = (J - 1) / HourlySteps
      Dim H = Limit(H_temp, 0, 23)
      NLDAS_RaPixels(J) = CalculateInstantaneousRa(
        NLDAS_Longitude,
        NLDAS_Latitude,
        NLDAS_Elevation,
        NLDAS_Slope,
        NLDAS_Aspect,
        NLDAS_AirTemperature,
        NLDAS_AirPressure,
        SolarPositions(J)
      )
    Next

    'Lapse Rate (F/ft)
    Dim M1 = RecordDate.Month - 2
    Dim M2 = RecordDate.Month - 1
    If RecordDate.Day >= 15 Then
      M1 = RecordDate.Month - 1
      M2 = RecordDate.Month
    End If
    If M1 = -1 Then M1 = 11
    If M2 = 12 Then M2 = 0
    Dim Period = DateTime.DaysInMonth(RecordDate.Year, M1 + 1)
    Dim Fraction = If(RecordDate.Day >= 15, RecordDate.Day - 15, Period - (15 - RecordDate.Day)) / Period
    Dim LapseRate = MonthlyLapseRate(M1) * (1 - Fraction) + MonthlyLapseRate(M2) * Fraction

    'Create Intermediate Arrays
    Dim FunctionDoY As Double = (RecordDate.DayOfYear - 1) / 365 * 2 * π
    Dim CosDoY As Double = Math.Cos(FunctionDoY)
    Dim SinDoY As Double = Math.Sin(FunctionDoY)

    'Sinusoidal Regression Correction Factors
    Dim TimeVariables(5) As Double
    TimeVariables(0) = 1
    TimeVariables(1) = CosDoY
    TimeVariables(2) = SinDoY

    Dim MinimumAirTemperature As Single = Single.MaxValue
    Dim MaximumAirTemperature As Single = Single.MinValue

    'For H = 0 To 23
    Dim H As Integer = 0

    'Record Date and Time
    Dim Time = RecordDate.AddHours(H)
    Dim FunctionHour As Double = (Time.Hour) / 23 * 2 * π
    TimeVariables(3) = Math.Cos(FunctionHour)
    TimeVariables(4) = Math.Sin(FunctionHour)

    'Correct for Topographical Differences
    Dim ΔZ = Elevation - NLDAS_Elevation

    'Air Temperature (F)
    AirTemperature += LapseRate * ΔZ

    TimeVariables(5) = AirTemperature
    AirTemperature -= SumProduct(AirTemperatureCorrectionCoefficients, TimeVariables)

    If AirTemperature < MinimumAirTemperature Then MinimumAirTemperature = AirTemperature
    If AirTemperature > MaximumAirTemperature Then MaximumAirTemperature = AirTemperature

    'Air Pressure (kPa)
    AirPressure = AdjustAirPressure(AirPressure, AirTemperature, ΔZ)
    'Return AirPressure

    'Solar Radiation (Langley/hr)
    Dim ProjectRa As Double = 0
    Dim Start As Integer = H * HourlySteps
    For J As Integer = Start To Start + HourlySteps
      ProjectRa += CalculateInstantaneousRa(
        Longitude, Latitude, Elevation, Slope, Aspect, AirTemperature, AirPressure, SolarPositions(J)
      )
    Next
    Dim HourlyAverage As Integer = HourlySteps + 1
    ProjectRa /= HourlyAverage

    Dim NLDAS_Ra As New Double
    For J As Integer = Start To Start + HourlySteps
      NLDAS_Ra += NLDAS_RaPixels(J)
    Next
    NLDAS_Ra /= HourlyAverage

    Dim SolarRadiation As Double = TranslateRs(NLDAS_Rs, NLDAS_Ra, ProjectRa, Slope)

    'Relative Humidity (%)
    Dim RelativeHumidity = NLDAS_RH
    TimeVariables(5) = RelativeHumidity
    RelativeHumidity -= SumProduct(HumidityCorrectionCoefficients, TimeVariables)
    If RelativeHumidity > 100 Then RelativeHumidity = 100
    If RelativeHumidity < 7 Then RelativeHumidity = 7

    'Wind Speed (mph)
    Dim WindSpeed = (NLDAS_WindU ^ 2 + NLDAS_WindV ^ 2) ^ 0.5
    If WindSpeed > 5.5 Then WindSpeed = 5.5

    'Calculate ASCE Reference Evapotranspiration (in)
    Dim ETo = CalculateHourlyASCEReferenceET(
      Elevation, AirTemperature, RelativeHumidity, WindSpeed,
      SolarRadiation, ProjectRa, ToFeet(2), ReferenceET.ShortReference, AirPressure
    )
    Return ETo

  End Function



  ''' <summary>
  ''' Calculates hourly reference evapotranspiration from the ASCE standardized equation.
  ''' </summary>
  ''' <param name="Elevation">Height Above Mean Sea Level (Feet)</param>
  ''' <param name="AirTemperature">Average Hourly Temperature (Fahrenheit)</param>
  ''' <param name="RelativeHumidity">Average Hourly Relative Humidity (Percent)</param>
  ''' <param name="SolarRadiation">Total Hourly Incoming Solar Radiation (Langleys)</param>
  ''' <param name="ExtraterrestrialRadiation">Total Hourly Extraterrestrial Solar Radiation (Langleys)</param>
  ''' <param name="WindSpeed">Average Hourly Wind Speed (Miles/Hour)</param>
  ''' <param name="AnemometerHeight">Height of Anemometer from the Ground (Feet)</param>
  ''' <param name="ReferenceET">Short(Grass) or Long(Alfalfa) Reference Height</param>
  ''' <param name="AirPressure">Average Hourly Air Pressure (kiloPascals)</param>
  ''' <returns>Estimated Reference Evapotranspiration (Inches/Hour)</returns>
  ''' <remarks>Source: ASCE Standardized Reference Evapotranspiration Equation (2005)</remarks>
  Function CalculateHourlyASCEReferenceET(ByVal Elevation As Double, ByVal AirTemperature As Double, ByVal RelativeHumidity As Double, ByVal WindSpeed As Double, ByVal SolarRadiation As Double, ByVal ExtraterrestrialRadiation As Double, ByVal AnemometerHeight As Double, ByVal ReferenceET As ReferenceET, Optional ByVal AirPressure As Double = Double.MinValue)
    'Mean Air Temperature (°C)
    Dim T As Double = ToCelsius(AirTemperature)

    'Elevation Above Mean Sea Level (m)
    Dim Z As Double = ToMeters(Elevation)

    'Atmospheric Pressure (kPa)
    Dim P As Double
    If AirPressure <> Double.MinValue Then
      P = AirPressure
    Else
      'P = 101.3 * ((293 - 0.0065 * Z) / 293) ^ 5.26
      P = (2.405994173 - 0.00005337529735 * Z) ^ 5.26
    End If

    'Psychrometric Constant (kPa/°C)
    Dim γ As Double = 0.000665 * P

    'Slope of the Saturation Vapor Pressure-Temperature Curve (kPa/°C)
    Dim Δ As Double = 2503 * Math.Exp(17.27 * T / (T + 237.3)) / (T + 237.3) ^ 2

    'Saturation Vapor Pressure (kPa)
    Dim Es As Double = 0.6108 * Math.Exp(17.27 * T / (T + 237.3))

    'Actual Vapor Pressure (kPa)
    Dim Ea As Double = 0.01 * RelativeHumidity * Es 'Conversion from Percent

    'Extraterrestrial Radiation (MJ/(m^2*h))
    Dim Ra As Double = ToMJperM2(ExtraterrestrialRadiation)

    'Calculated Clear-sky Radiation (MJ/(m^2*h))
    Dim Rso As Double = (0.75 + 0.00002 * Z) * Ra

    'Measured Solar Radiation (MJ/(m^2*h))
    Dim Rs As Double = ToMJperM2(SolarRadiation)

    'Albedo (unitless)
    Dim α As Double = 0.23

    'Net Outgoing Long-Wave Radiation (MJ/(m^2*h))
    Dim Rns As Double = (1 - α) * Rs

    'Cloudiness Function (unitless)
    Dim Fcd As Double = 0.7
    If Rso > 0 Then
      Dim Temp As Double = Rs / Rso
      Fcd = 1.35 * Limit(Temp, 0.3, 1) - 0.35
    End If

    'Stefan-Boltzmann Constant (MJ/(K^4*m^2*h))
    Dim σ As Double = 0.0000000002042

    'Net Shortwave Radiation (MJ/(m^2*h))
    Dim Rnl As Double = σ * Fcd * (0.34 - 0.14 * Ea ^ 0.5) * ToKelvin(T) ^ 4

    'Net Radiation (MJ/(m^2*h))
    Dim Rn As Double = Rns - Rnl

    'Soil Heat Flux Density (MJ/(m^2*h))
    Dim G As Double

    'Numerator Constant that Changes with Reference Type and Calculation Time Step (K*mm*s^3/(Mg*h))
    Dim Cn As Double

    'Denominator Constant that Changes with Reference Type and Calculation Time Step (s/m)
    Dim Cd As Double

    'Assignment of Reference Values
    If Rn >= 0 Then
      Select Case ReferenceET
        Case ReferenceET.ShortReference
          G = 0.1 * Rn
          Cn = 37
          Cd = 0.24
        Case ReferenceET.LongReference
          G = 0.04 * Rn
          Cn = 66
          Cd = 0.25
      End Select
    Else
      Select Case ReferenceET
        Case ReferenceET.ShortReference
          G = 0.5 * Rn
          Cn = 37
          Cd = 0.96
        Case ReferenceET.LongReference
          G = 0.2 * Rn
          Cn = 66
          Cd = 1.7
      End Select
    End If

    'Height of Wind Measurement Above Ground Surface (m)
    Dim Zw As Double = ToMeters(AnemometerHeight)

    'Measured Wind Speed at Zw (m/s)
    Dim Uz As Double = ToMetersPerSecond(WindSpeed)

    'Adjusted Wind Speed at 2 m Above Ground Surface (m/s)
    Dim U2 As Double = 4.87 * Uz / Math.Log(67.8 * Zw - 5.42)

    'Inverse Latent Heat of Vaporization (kg/MJ)
    Dim λ As Double = 0.408

    'ASCE Standardized Reference Evapotranspiration (mm/hr)
    Dim ET As Double = (λ * Δ * (Rn - G) + γ * Cn * U2 * (Es - Ea) / (T + 273)) / (Δ + γ * (1 + Cd * U2))

    Return ToInches(ET)
  End Function


  ''' <summary>
  ''' Calculates instantaneous incoming extraterrestrial radiation for a given location and time period on the earth for a sloped surface.
  ''' </summary>
  ''' <param name="Longitude">(Decimal Degrees)</param>
  ''' <param name="Latitude">(Decimal Degrees)</param>
  ''' <param name="Elevation">Height Above Mean Sea Level (Feet)</param>
  ''' <param name="Slope">Positive Inclination from Horizontal (Decimal Degrees)</param>
  ''' <param name="Azimuth">Positive Clockwise Slope Orientation from North(0) [East(90), South(180), West(270)] (Decimal Degrees)</param>
  ''' <param name="Temperature">Air Temperature (Fahrenheit)</param>
  ''' <param name="Pressure">Air Pressure (kPa)</param>
  ''' <param name="SolarPosition">Earth-Sun Orientation and Distance</param>
  ''' <returns>Estimated Extraterrestrial Radiation for a Sloped Surface (Langleys/Hour)</returns>
  Function CalculateInstantaneousRa(ByVal Longitude As Double, ByVal Latitude As Double, ByVal Elevation As Double, ByVal Slope As Double, ByVal Azimuth As Double, ByVal Temperature As Double, ByVal Pressure As Double, ByVal SolarPosition As SolarPosition) As Double
    'Variable Conversions
    Dim σ As Double = ToRadians(Longitude)
    Dim ϕ As Double = ToRadians(Latitude)
    Dim Z As Double = ToMeters(Elevation)
    Dim ω As Double = ToRadians(Slope)
    Dim γ As Double = ToRadians(Azimuth - 180)

    Dim Cosϕ As Double = Math.Cos(ϕ)
    Dim Sinϕ As Double = Math.Sin(ϕ)

    'Observer Local Hour Angle (radians)
    Dim H As Double = LimitAngle(SolarPosition.ν + σ - SolarPosition.α)
    Dim CosH As Double = Math.Cos(H)

    'Equatorial Horizontal Parallax of the Sun (radians)
    Dim ξ As Double = ToRadians(8.794 / (3600 * SolarPosition.R))
    Dim Sinξ As Double = Math.Sin(ξ)

    'Parallax in the Sun Right Ascension
    Dim U As Double = Math.Atan(0.99664719 * Math.Tan(ϕ))
    Dim X As Double = Math.Cos(U) + Z * Cosϕ / 6378140
    Dim Y As Double = 0.99664719 * Math.Sin(U) + Z * Sinϕ / 6378140
    Dim Δα As Double = Math.Atan2(-X * Sinξ * Math.Sin(H), Math.Cos(SolarPosition.δ) - X * Sinξ * CosH)

    'Topocentric Sun Right Ascension (radians)
    'Dim αp As Double = SolarPosition.α + Δα

    'Topocentric Sun Declination (radians)
    Dim δp As Double = Math.Atan2((Math.Sin(SolarPosition.δ) - Y * Sinξ) * Math.Cos(Δα), Math.Cos(SolarPosition.δ) - X * Sinξ * CosH)

    'Topocentric Local Hour Angle (radians)
    Dim Hp As Double = H - Δα
    Dim CosHp As Double = Math.Cos(Hp)

    'Topocentric Elevation Angle Without Atmospheric Refraction (radians)
    Dim e0 As Double = Math.Asin(Sinϕ * Math.Sin(δp) + Cosϕ * Math.Cos(δp) * CosHp)

    'Atmospheric Refraction Correction (radians)
    Dim Δe As Double = ToRadians((Pressure / 101) * (283 / (273 + ToCelsius(Temperature))) * 1.02 / (60 * Math.Tan(ToRadians(ToDegrees(e0) + 10.3 / (ToDegrees(e0) + 5.11)))))

    'Topocentric Elevation Angle (radians)
    Dim e As Double = e0 + Δe

    'Zenith Angle (radians)
    Dim θ As Double = π / 2 - e

    'Topocentric Astronomers Azimuth Angle (radians)
    Dim Γc As Double = Math.Atan2(Math.Sin(Hp), CosHp * Sinϕ - Math.Tan(δp) * Cosϕ)

    'Topocentric Azimuth Angle (radians)
    Dim Φ As Double = LimitAngle(Γc + π)

    'Incidence Angle for a Surface Oriented in Any Direction
    Dim I As Double = Math.Acos(Math.Cos(θ) * Math.Cos(ω) + Math.Sin(ω) * Math.Sin(θ) * Math.Cos(Γc - γ))

    'Solar Constant (W/m^2)
    Dim Gsc As Double = 1367

    'Extraterrestrial Radiation (W/m^2)
    Dim Ra As Double = 0
    If e > 0 And π / 2 - I > 0 Then Ra = Gsc * Math.Cos(I) / SolarPosition.R ^ 2

    Return ToLangleysPerHour(Ra)
  End Function


  ''' <summary>
  ''' Translates incoming solar radiation from one site to another taking into account slope differences (Source interpolation site assumed nearly horizontal).
  ''' </summary>
  ''' <param name="RsSource">Incoming Solar Radiation at Source Interpolation Site (Same Units as Other Input Radiations)</param>
  ''' <param name="RaSource">Extraterrestrial Solar Radiation at Source Interpolation Site (Same Units as Other Input Radiations)</param>
  ''' <param name="RaDestination">Extraterrestrial Solar Radiation at Source Interpolation Site (Same Units as Other Input Radiations)</param>
  ''' <param name="Slope">Positive Inclination from Horizontal (Decimal Degrees)</param>
  ''' <returns>Incoming Solar Radiation at Destination Translation Site (Same Units as Input Radiations)</returns>
  ''' <remarks>Adapted from: Allen, R. G., Trezza, R., Tasumi, M. (2006). Analytical integrated functions for daily solar radiation on slopes. Agricultural and Forest Meteorology, 139(1), 55-73.</remarks>
  Function TranslateRs(ByVal RsSource As Double, ByVal RaSource As Double, ByVal RaDestination As Double, ByVal Slope As Double) As Double
    If RsSource <= 0 Or RaSource <= 0 Then
      Return 0
    Else
      'Convert to Radians
      Slope = ToRadians(Slope)

      'Actual Atmospheric Transmissivity (Unitless)
      Dim Tau_sw As Double = RsSource / RaSource

      'Source Actual Direct Radiation Index (Unitless)
      Dim Kb As Double = 0
      Select Case Tau_sw
        Case Is >= 0.42
          Kb = 1.56 * Tau_sw - 0.55
        Case Is <= 0.175
          Kb = 0.016 * Tau_sw
        Case Else
          Kb = ((0.765 * Tau_sw + 0.828) * Tau_sw + 0.28) * Tau_sw + 0.022
      End Select

      'Source Actual Diffuse Radiation Index (Unitless)
      Dim Kd As Double = Tau_sw - Kb

      'Direct Radiation Fraction (Unitless)
      Dim Fb As Double = RaDestination / RaSource

      'Reflectance Factor (Unitless)
      Dim Fi As Double = 0.75 + 0.25 * Math.Cos(Slope) - 0.5 * Slope / π

      'Sky-view Factor (Unitless)
      Dim Fia As Double = (1 - Kb) * (1 + (Kb / (Kb + Kd)) ^ 0.5 * Math.Sin(Slope / 2) ^ 3) * Fi + Fb * Kb

      'Incoming Solar Radiation at Interpolated Site (Source Site Units)
      Dim RsDestination As Double = RsSource * ((Fb * Kb + Fia * Kd) / tau_sw + 0.23 * (1 - Fi))

      Return RsDestination
    End If
  End Function


  ''' <summary>
  ''' Calculates the geocentric sun angles and distance at a given instant in time.
  ''' </summary>
  ''' <param name="RecordDate">Event Time</param>
  ''' <returns>Earth-Sun Orientation and Distance</returns>
  Function CalculateSolarPosition(ByVal RecordDate As DateTime) As SolarPosition

    'Julian Day
    Dim Month As Integer = RecordDate.Month
    Dim Year As Integer = RecordDate.Year
    If Month < 3 Then
      Month += 12
      Year -= 1
    End If
    Dim JD As Double = Int(365.25 * (Year + 4716)) + Int(30.6001 * (Month + 1)) + RecordDate.Day + RecordDate.TimeOfDay.TotalHours / 24 - 1524.5
    If JD > 2299160 Then
      Dim A = Int(Year / 100)
      JD += (2 - A + Int(A / 4))
    End If

    ' Approximate Difference between the Earth Rotation Time and Terrestrial Time (seconds)
    Dim ΔT As Double = LookupΔT(RecordDate)

    'Julian Day, Ephemeris
    Dim JDE As Double = JD + ΔT / 86400 'seconds to day fraction

    'Julian Century, Ephemeris
    Dim JC As Double = (JD - 2451545) / 36525
    Dim JCE As Double = (JDE - 2451545) / 36525

    'Julian Millennium, Ephemeris
    Dim JME As Double = JCE / 10

    'Earth Heliocentric Longitude (radians)
    Dim L As Double = CalculateEarthHeliocentricProperty(JME, EarthHeliocentricType.Longitude)

    'Earth Heliocentric Latitude (radians)
    Dim B As Double = CalculateEarthHeliocentricProperty(JME, EarthHeliocentricType.Latitude)

    'Earth-Sun Distance (astronomical units, au)
    Dim R As Double = CalculateEarthHeliocentricProperty(JME, EarthHeliocentricType.Radius)

    'Geocentric Longitude (radians)
    Dim Θ As Double = LimitAngle(L + π)

    'Geocentric Latitude (radians)
    Dim β As Double = -B

    'Mean Elongation of the Moon from the Sun (degrees)
    Dim X(4) As Double
    X(0) = CalculateSolarPolynomial(JCE, 1 / 189474, -0.0019142, 445267.11148, 297.85036)

    'Mean Anomaly of the Sun-Earth (degrees)
    X(1) = CalculateSolarPolynomial(JCE, 1 / -300000, -0.0001603, 35999.05034, 357.52772)

    'Mean Anomaly of the Moon (degrees)
    X(2) = CalculateSolarPolynomial(JCE, 1 / 56250, 0.0086972, 477198.867398, 134.96298)

    'Moon's Argument of Latitude  (degrees)
    X(3) = CalculateSolarPolynomial(JCE, 1 / 327270, -0.0036825, 483202.017538, 93.27191)

    'Longitude of the Ascending Node of the Moon’s Mean Orbit on the Ecliptic, Measured from the Mean Equinox of the Date (degrees)
    X(4) = CalculateSolarPolynomial(JCE, 1 / 450000, 0.0020708, -1934.136261, 125.04452)

    'Nutation in Longitude (radians)
    Dim Δψ As Double = CalculateNutation(JCE, X, NutationType.Longitude)

    'Nutation In Obliquity (radians)
    Dim Δε As Double = CalculateNutation(JCE, X, NutationType.Obliquity)

    'Mean Obliquity of the Ecliptic (arc seconds)
    Dim U As Double = JME / 10
    Dim ε0 As Double = 84381.448 + U * (-4680.96 + U * (-1.55 + U * (1999.25 + U * (-51.38 + U * (-249.67 + U * (-39.05 + U * (7.12 + U * (27.87 + U * (5.79 + U * 2.45)))))))))

    'True Obliquity of the Ecliptic (radians)
    Dim ε As Double = ToRadians(ε0 / 3600) + Δε

    'Aberration Correction (radians)
    Dim Δτ As Double = ToRadians(-20.4898 / (3600 * R))

    'Apparent Sun Longitude (radians)
    Dim λ As Double = Θ + Δψ + Δτ

    'Mean Sidereal Time at Greenwich (radians)
    Dim ν0 As Double = LimitAngle(ToRadians(280.46061837 + 360.98564736629 * (JD - 2451545) + JC * JC * (0.000387933 - JC / 38710000)))

    'Trigonometric Values
    Dim Cosε As Double = Math.Cos(ε)
    Dim Sinε As Double = Math.Sin(ε)
    Dim Sinλ As Double = Math.Sin(λ)

    'Apparent Sidereal Time at Greenwich (radians)
    Dim ν As Double = ν0 + Δψ * Cosε

    'Geocentric Sun Right Ascension (radians)
    Dim α As Double = LimitAngle(Math.Atan2(Sinλ * Cosε - Math.Tan(β) * Sinε, Math.Cos(λ)))

    'Geocentric Sun Declination (radians)
    Dim δ As Double = Math.Asin(Math.Sin(β) * Cosε + Math.Cos(β) * Sinε * Sinλ)

    Return New SolarPosition(R, δ, α, ν)
  End Function


  ''' <summary>
  ''' Returns a yearly average difference between universal time and actual terrestrial time (1620-2014).
  ''' </summary>
  ''' <param name="RecordDate">Event Time</param>
  ''' <returns>Time Difference (seconds)</returns>
  Function LookupΔT(ByVal RecordDate As DateTime) As Double
    Dim ΔT = New Double() {124, 119, 115, 110, 106, 102, 98, 95, 91, 88, 85, 82, 79, 77, 74, 72, 70, 67, 65, 63, 62, 60, 58, 57, 55, 54, 53, 51, 50, 49, 48, 47, 46, 45, 44, 43, 42, 41, 40, 38, 37, 36, 35, 34, 33, 32, 31, 30, 28, 27, 26, 25, 24, 23, 22, 21, 20, 19, 18, 17, 16, 15, 14, 14, 13, 12, 12, 11, 11, 10, 10, 10, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 10, 10, 10, 10, 10, 10, 10, 10, 10, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 13, 13, 13, 13, 13, 13, 13, 14, 14, 14, 14, 14, 14, 14, 15, 15, 15, 15, 15, 15, 15, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 16, 16, 16, 16, 15, 15, 14, 14, 13.7, 13.4, 13.1, 12.9, 12.7, 12.6, 12.5, 12.5, 12.5, 12.5, 12.5, 12.5, 12.5, 12.5, 12.5, 12.5, 12.5, 12.4, 12.3, 12.2, 12, 11.7, 11.4, 11.1, 10.6, 10.2, 9.6, 9.1, 8.6, 8, 7.5, 7, 6.6, 6.3, 6, 5.8, 5.7, 5.6, 5.6, 5.6, 5.7, 5.8, 5.9, 6.1, 6.2, 6.3, 6.5, 6.6, 6.8, 6.9, 7.1, 7.2, 7.3, 7.4, 7.5, 7.6, 7.7, 7.7, 7.8, 7.8, 7.88, 7.82, 7.54, 6.97, 6.4, 6.02, 5.41, 4.1, 2.92, 1.82, 1.61, 0.1, -1.02, -1.28, -2.69, -3.24, -3.64, -4.54, -4.71, -5.11, -5.4, -5.42, -5.2, -5.46, -5.46, -5.79, -5.63, -5.64, -5.8, -5.66, -5.87, -6.01, -6.19, -6.64, -6.44, -6.47, -6.09, -5.76, -4.66, -3.74, -2.72, -1.54, -0.02, 1.24, 2.64, 3.86, 5.37, 6.14, 7.75, 9.13, 10.46, 11.53, 13.36, 14.65, 16.01, 17.2, 18.24, 19.06, 20.25, 20.95, 21.16, 22.25, 22.41, 23.03, 23.49, 23.62, 23.86, 24.49, 24.34, 24.08, 24.02, 24, 23.87, 23.95, 23.86, 23.93, 23.73, 23.92, 23.96, 24.02, 24.33, 24.83, 25.3, 25.7, 26.24, 26.77, 27.28, 27.78, 28.25, 28.71, 29.15, 29.57, 29.97, 30.36, 30.72, 31.07, 31.35, 31.68, 32.18, 32.68, 33.15, 33.59, 34, 34.47, 35.03, 35.73, 36.54, 37.43, 38.29, 39.2, 40.18, 41.17, 42.23, 43.37, 44.49, 45.48, 46.46, 47.52, 48.53, 49.59, 50.54, 51.38, 52.17, 52.96, 53.79, 54.34, 54.87, 55.32, 55.82, 56.3, 56.86, 57.57, 58.31, 59.12, 59.99, 60.78, 61.63, 62.3, 62.97, 63.47, 63.83, 64.09, 64.3, 64.47, 64.57, 64.69, 64.85, 65.15, 65.46, 65.78, 66.2, 66.45, 66.74, 67.09, 67.43}

    Dim Y_Index As Double = RecordDate.Year - 1620
    Dim Y_Min As Double = 0
    Dim Y_Max As Double = ΔT.Length - 1
    Dim Y = Limit(Y_Index, Y_Min, Y_Max)

    Return ΔT(Y)
  End Function


  ''' <summary>
  ''' Calculates earth heliocentric longitude, latitude, or earth-sun distance.
  ''' </summary>
  ''' <param name="JM">Julian Millenium</param>
  ''' <param name="Type">Earth Heliocentric Property Type</param>
  ''' <returns>Earth Heliocentric Property [Longitude or Latitude (Decimal Degrees), or Earth-Sun Distance (astronomical units)]</returns>
  Function CalculateEarthHeliocentricProperty(ByVal JM As Double, ByVal Type As EarthHeliocentricType) As Double
    Dim ABC()(,) As Double = Nothing
    Select Case Type
      Case EarthHeliocentricType.Longitude
        ABC = New Double()(,) {New Double(,) {{175347046, 0, 0}, {3341656, 4.6692568, 6283.07585}, {34894, 4.6261, 12566.1517}, {3497, 2.7441, 5753.3849}, {3418, 2.8289, 3.5231}, {3136, 3.6277, 77713.7715}, {2676, 4.4181, 7860.4194}, {2343, 6.1352, 3930.2097}, {1324, 0.7425, 11506.7698}, {1273, 2.0371, 529.691}, {1199, 1.1096, 1577.3435}, {990, 5.233, 5884.927}, {902, 2.045, 26.298}, {857, 3.508, 398.149}, {780, 1.179, 5223.694}, {753, 2.533, 5507.553}, {505, 4.583, 18849.228}, {492, 4.205, 775.523}, {357, 2.92, 0.067}, {317, 5.849, 11790.629}, {284, 1.899, 796.298}, {271, 0.315, 10977.079}, {243, 0.345, 5486.778}, {206, 4.806, 2544.314}, {205, 1.869, 5573.143}, {202, 2.458, 6069.777}, {156, 0.833, 213.299}, {132, 3.411, 2942.463}, {126, 1.083, 20.775}, {115, 0.645, 0.98}, {103, 0.636, 4694.003}, {102, 0.976, 15720.839}, {102, 4.267, 7.114}, {99, 6.21, 2146.17}, {98, 0.68, 155.42}, {86, 5.98, 161000.69}, {85, 1.3, 6275.96}, {85, 3.67, 71430.7}, {80, 1.81, 17260.15}, {79, 3.04, 12036.46}, {75, 1.76, 5088.63}, {74, 3.5, 3154.69}, {74, 4.68, 801.82}, {70, 0.83, 9437.76}, {62, 3.98, 8827.39}, {61, 1.82, 7084.9}, {57, 2.78, 6286.6}, {56, 4.39, 14143.5}, {56, 3.47, 6279.55}, {52, 0.19, 12139.55}, {52, 1.33, 1748.02}, {51, 0.28, 5856.48}, {49, 0.49, 1194.45}, {41, 5.37, 8429.24}, {41, 2.4, 19651.05}, {39, 6.17, 10447.39}, {37, 6.04, 10213.29}, {37, 2.57, 1059.38}, {36, 1.71, 2352.87}, {36, 1.78, 6812.77}, {33, 0.59, 17789.85}, {30, 0.44, 83996.85}, {30, 2.74, 1349.87}, {25, 3.16, 4690.48}}, New Double(,) {{628331966747, 0, 0}, {206059, 2.678235, 6283.07585}, {4303, 2.6351, 12566.1517}, {425, 1.59, 3.523}, {119, 5.796, 26.298}, {109, 2.966, 1577.344}, {93, 2.59, 18849.23}, {72, 1.14, 529.69}, {68, 1.87, 398.15}, {67, 4.41, 5507.55}, {59, 2.89, 5223.69}, {56, 2.17, 155.42}, {45, 0.4, 796.3}, {36, 0.47, 775.52}, {29, 2.65, 7.11}, {21, 5.34, 0.98}, {19, 1.85, 5486.78}, {19, 4.97, 213.3}, {17, 2.99, 6275.96}, {16, 0.03, 2544.31}, {16, 1.43, 2146.17}, {15, 1.21, 10977.08}, {12, 2.83, 1748.02}, {12, 3.26, 5088.63}, {12, 5.27, 1194.45}, {12, 2.08, 4694}, {11, 0.77, 553.57}, {10, 1.3, 6286.6}, {10, 4.24, 1349.87}, {9, 2.7, 242.73}, {9, 5.64, 951.72}, {8, 5.3, 2352.87}, {6, 2.65, 9437.76}, {6, 4.67, 4690.48}}, New Double(,) {{52919, 0, 0}, {8720, 1.0721, 6283.0758}, {309, 0.867, 12566.152}, {27, 0.05, 3.52}, {16, 5.19, 26.3}, {16, 3.68, 155.42}, {10, 0.76, 18849.23}, {9, 2.06, 77713.77}, {7, 0.83, 775.52}, {5, 4.66, 1577.34}, {4, 1.03, 7.11}, {4, 3.44, 5573.14}, {3, 5.14, 796.3}, {3, 6.05, 5507.55}, {3, 1.19, 242.73}, {3, 6.12, 529.69}, {3, 0.31, 398.15}, {3, 2.28, 553.57}, {2, 4.38, 5223.69}, {2, 3.75, 0.98}}, New Double(,) {{289, 5.844, 6283.076}, {35, 0, 0}, {17, 5.49, 12566.15}, {3, 5.2, 155.42}, {1, 4.72, 3.52}, {1, 5.3, 18849.23}, {1, 5.97, 242.73}}, New Double(,) {{114, 3.142, 0}, {8, 4.13, 6283.08}, {1, 3.84, 12566.15}}, New Double(,) {{1, 3.14, 0}}}
      Case EarthHeliocentricType.Latitude
        ABC = New Double()(,) {New Double(,) {{280, 3.199, 84334.662}, {102, 5.422, 5507.553}, {80, 3.88, 5223.69}, {44, 3.7, 2352.87}, {32, 4, 1577.34}}, New Double(,) {{9, 3.9, 5507.55}, {6, 1.73, 5223.69}}}
      Case EarthHeliocentricType.Radius
        ABC = New Double()(,) {New Double(,) {{100013989, 0, 0}, {1670700, 3.0984635, 6283.07585}, {13956, 3.05525, 12566.1517}, {3084, 5.1985, 77713.7715}, {1628, 1.1739, 5753.3849}, {1576, 2.8469, 7860.4194}, {925, 5.453, 11506.77}, {542, 4.564, 3930.21}, {472, 3.661, 5884.927}, {346, 0.964, 5507.553}, {329, 5.9, 5223.694}, {307, 0.299, 5573.143}, {243, 4.273, 11790.629}, {212, 5.847, 1577.344}, {186, 5.022, 10977.079}, {175, 3.012, 18849.228}, {110, 5.055, 5486.778}, {98, 0.89, 6069.78}, {86, 5.69, 15720.84}, {86, 1.27, 161000.69}, {65, 0.27, 17260.15}, {63, 0.92, 529.69}, {57, 2.01, 83996.85}, {56, 5.24, 71430.7}, {49, 3.25, 2544.31}, {47, 2.58, 775.52}, {45, 5.54, 9437.76}, {43, 6.01, 6275.96}, {39, 5.36, 4694}, {38, 2.39, 8827.39}, {37, 0.83, 19651.05}, {37, 4.9, 12139.55}, {36, 1.67, 12036.46}, {35, 1.84, 2942.46}, {33, 0.24, 7084.9}, {32, 0.18, 5088.63}, {32, 1.78, 398.15}, {28, 1.21, 6286.6}, {28, 1.9, 6279.55}, {26, 4.59, 10447.39}}, New Double(,) {{103019, 1.10749, 6283.07585}, {1721, 1.0644, 12566.1517}, {702, 3.142, 0}, {32, 1.02, 18849.23}, {31, 2.84, 5507.55}, {25, 1.32, 5223.69}, {18, 1.42, 1577.34}, {10, 5.91, 10977.08}, {9, 1.42, 6275.96}, {9, 0.27, 5486.78}}, New Double(,) {{4359, 5.7846, 6283.0758}, {124, 5.579, 12566.152}, {12, 3.14, 0}, {9, 3.63, 77713.77}, {6, 1.87, 5573.14}, {3, 5.47, 18849.23}}, New Double(,) {{145, 4.273, 6283.076}, {7, 3.92, 12566.15}}, New Double(,) {{4, 2.56, 6283.08}}}
    End Select

    Dim Value(ABC.GetLength(0) - 1) As Double
    For I As Integer = 0 To ABC.GetLength(0) - 1
      For J As Integer = 0 To ABC(I).GetLength(0) - 1
        Value(I) += ABC(I)(J, 0) * Math.Cos(ABC(I)(J, 1) + ABC(I)(J, 2) * JM)
      Next
    Next

    Dim EarthHeliocentricValue As Double = 0
    Dim JMpower As Double = 1
    For I As Integer = 0 To Value.Length - 1
      EarthHeliocentricValue += Value(I) * JMpower
      JMpower *= JM
    Next
    EarthHeliocentricValue = EarthHeliocentricValue / 10 ^ 8
    If Type = EarthHeliocentricType.Longitude Then EarthHeliocentricValue = LimitAngle(EarthHeliocentricValue)

    Return EarthHeliocentricValue
  End Function


  ''' <summary>
  ''' Calculates a third order polynomial for earth solar calculations (A*V^3+B*V^2+C*V+D).
  ''' </summary>
  ''' <param name="Variable">Input Variable Raised to a Power</param>
  ''' <param name="A">3rd Order Constant</param>
  ''' <param name="B">2nd Order Constant</param>
  ''' <param name="C">1st Order Constant</param>
  ''' <param name="D">Constant</param>
  ''' <returns>Calculated Value</returns>
  Function CalculateSolarPolynomial(ByVal Variable As Double, ByVal A As Double, ByVal B As Double, ByVal C As Double, ByVal D As Double) As Double
    Return LimitAngle(((A * Variable + B) * Variable + C) * Variable + D, 360)
  End Function


  ''' <summary>
  ''' Calculates nutation in longitude or obliquity for an input Julian Century.
  ''' </summary>
  ''' <param name="JC">Julian Century</param>
  ''' <param name="X">Earth-Moon-Sun Properties Array</param>
  ''' <param name="Type">Nutation Type</param>
  ''' <returns>Nutation</returns>
  Function CalculateNutation(ByVal JC As Double, ByVal X() As Double, ByVal Type As NutationType) As Double
    Dim Y = New Double(,) {{0, 0, 0, 0, 1}, {-2, 0, 0, 2, 2}, {0, 0, 0, 2, 2}, {0, 0, 0, 0, 2}, {0, 1, 0, 0, 0}, {0, 0, 1, 0, 0}, {-2, 1, 0, 2, 2}, {0, 0, 0, 2, 1}, {0, 0, 1, 2, 2}, {-2, -1, 0, 2, 2}, {-2, 0, 1, 0, 0}, {-2, 0, 0, 2, 1}, {0, 0, -1, 2, 2}, {2, 0, 0, 0, 0}, {0, 0, 1, 0, 1}, {2, 0, -1, 2, 2}, {0, 0, -1, 0, 1}, {0, 0, 1, 2, 1}, {-2, 0, 2, 0, 0}, {0, 0, -2, 2, 1}, {2, 0, 0, 2, 2}, {0, 0, 2, 2, 2}, {0, 0, 2, 0, 0}, {-2, 0, 1, 2, 2}, {0, 0, 0, 2, 0}, {-2, 0, 0, 2, 0}, {0, 0, -1, 2, 1}, {0, 2, 0, 0, 0}, {2, 0, -1, 0, 1}, {-2, 2, 0, 2, 2}, {0, 1, 0, 0, 1}, {-2, 0, 1, 0, 1}, {0, -1, 0, 0, 1}, {0, 0, 2, -2, 0}, {2, 0, -1, 2, 1}, {2, 0, 1, 2, 2}, {0, 1, 0, 2, 2}, {-2, 1, 1, 0, 0}, {0, -1, 0, 2, 2}, {2, 0, 0, 2, 1}, {2, 0, 1, 0, 0}, {-2, 0, 2, 2, 2}, {-2, 0, 1, 2, 1}, {2, 0, -2, 0, 1}, {2, 0, 0, 0, 1}, {0, -1, 1, 0, 0}, {-2, -1, 0, 2, 1}, {-2, 0, 0, 0, 1}, {0, 0, 2, 2, 1}, {-2, 0, 2, 0, 1}, {-2, 1, 0, 2, 1}, {0, 0, 1, -2, 0}, {-1, 0, 1, 0, 0}, {-2, 1, 0, 0, 0}, {1, 0, 0, 0, 0}, {0, 0, 1, 2, 0}, {0, 0, -2, 2, 2}, {-1, -1, 1, 0, 0}, {0, 1, 1, 0, 0}, {0, -1, 1, 2, 2}, {2, -1, -1, 2, 2}, {0, 0, 3, 2, 2}, {2, -1, 0, 2, 2}}
    Dim ABCD = New Double(,) {{-171996, -174.2, 92025, 8.9}, {-13187, -1.6, 5736, -3.1}, {-2274, -0.2, 977, -0.5}, {2062, 0.2, -895, 0.5}, {1426, -3.4, 54, -0.1}, {712, 0.1, -7, 0}, {-517, 1.2, 224, -0.6}, {-386, -0.4, 200, 0}, {-301, 0, 129, -0.1}, {217, -0.5, -95, 0.3}, {-158, 0, 0, 0}, {129, 0.1, -70, 0}, {123, 0, -53, 0}, {63, 0, 0, 0}, {63, 0.1, -33, 0}, {-59, 0, 26, 0}, {-58, -0.1, 32, 0}, {-51, 0, 27, 0}, {48, 0, 0, 0}, {46, 0, -24, 0}, {-38, 0, 16, 0}, {-31, 0, 13, 0}, {29, 0, 0, 0}, {29, 0, -12, 0}, {26, 0, 0, 0}, {-22, 0, 0, 0}, {21, 0, -10, 0}, {17, -0.1, 0, 0}, {16, 0, -8, 0}, {-16, 0.1, 7, 0}, {-15, 0, 9, 0}, {-13, 0, 7, 0}, {-12, 0, 6, 0}, {11, 0, 0, 0}, {-10, 0, 5, 0}, {-8, 0, 3, 0}, {7, 0, -3, 0}, {-7, 0, 0, 0}, {-7, 0, 3, 0}, {-7, 0, 3, 0}, {6, 0, 0, 0}, {6, 0, -3, 0}, {6, 0, -3, 0}, {-6, 0, 3, 0}, {-6, 0, 3, 0}, {5, 0, 0, 0}, {-5, 0, 3, 0}, {-5, 0, 3, 0}, {-5, 0, 3, 0}, {4, 0, 0, 0}, {4, 0, 0, 0}, {4, 0, 0, 0}, {-4, 0, 0, 0}, {-4, 0, 0, 0}, {-4, 0, 0, 0}, {3, 0, 0, 0}, {-3, 0, 0, 0}, {-3, 0, 0, 0}, {-3, 0, 0, 0}, {-3, 0, 0, 0}, {-3, 0, 0, 0}, {-3, 0, 0, 0}, {-3, 0, 0, 0}}

    Select Case Type
      Case NutationType.Longitude
        Dim Δψ As Double = 0
        For I As Integer = 0 To Y.GetLength(0) - 1
          Dim ΣXY As Double = 0
          For J As Integer = 0 To Y.GetLength(1) - 1
            ΣXY += X(J) * Y(I, J)
          Next
          Δψ += (ABCD(I, 0) + ABCD(I, 1) * JC) * Math.Sin(ToRadians(ΣXY))
        Next
        CalculateNutation = Δψ
      Case NutationType.Obliquity
        Dim Δε As Double = 0
        For I As Integer = 0 To Y.GetLength(0) - 1
          Dim ΣXY As Double = 0
          For J As Integer = 0 To Y.GetLength(1) - 1
            ΣXY += X(J) * Y(I, J)
          Next
          Δε += (ABCD(I, 2) + ABCD(I, 3) * JC) * Math.Cos(ToRadians(ΣXY))
        Next
        CalculateNutation = Δε
      Case Else
        CalculateNutation = Double.MinValue
    End Select

    Return ToRadians(CalculateNutation / 36000000)
  End Function


  ''' <summary>
  ''' Redefines an angle between 0 to 2π or a specified upper boundary.
  ''' </summary>
  ''' <param name="Value">Input angle (Radians [or User-Defined])</param>
  ''' <param name="UpperBounds">Default = 2π [or User-Defined]</param>
  ''' <returns>Basic Positive Angle</returns>
  Function LimitAngle(ByVal Value As Double, Optional ByVal UpperBounds As Double = 2 * π) As Double
    Dim Output As Double = Value Mod UpperBounds

    If Output < 0 Then Output += UpperBounds

    Return Output
  End Function


  Function SumProduct(ByVal Array1 As Double(), ByVal Array2() As Double) As Double
    Dim Sum As Double = 0

    For I As Integer = 0 To Array1.Length - 1
        Sum += Array1(I) * Array2(I)
    Next

    Return Sum
  End Function


  ''' <summary>
  ''' Adjusts wind speed from one height to another using a logarithmic profile.
  ''' </summary>
  ''' <param name="Zw">Wind reference height (ft)</param>
  ''' <param name="H">Vegetation height (ft)</param>
  ''' <param name="Zu">Desired wind reference height (ft)</param>
  ''' <returns>Factor to adjust wind speed</returns>
  ''' <remarks></remarks>
  Function CalculateWindSpeedAdjustmentFactor(ByVal Zw As Double, Optional ByVal H As Double = 0.12 / 0.3048, Optional ByVal Zu As Double = 2 / 0.3048)
    Dim D = 0.67 * H 'zero plane displacement height
    Dim Zom = 0.123 * H 'aerodynamic roughness length
    Return Math.Log((Zu - D) / Zom) / Math.Log((Zw - D) / Zom) 'adjustment factor
  End Function


  ''' <summary>
  ''' Adjusts pressure for change in elevation.
  ''' </summary>
  ''' <param name="AirPressure">Average Air Pressure at Referenced Site (Any Unit--Multiplier Relative)</param>
  ''' <param name="AirTemperature">Average Air Temperature at Interpolated Site (Fahrenheit)</param>
  ''' <param name="ΔZ">Difference in Elevation of Referenced Site from Interpolated Site (Feet)</param>
  ''' <returns></returns>
  ''' <remarks></remarks>
  Function AdjustAirPressure(ByVal AirPressure As Double, ByVal AirTemperature As Double, ByVal ΔZ As Double) As Double
    'Acceleration of Gravity on Earth (ft/s^2)
    Dim g As Double = 32.174

    'Gas Constant for Air (ft-lb/(slug-R))
    Dim R As Double = 1716

    Return AirPressure / Math.Exp(g * ΔZ / (R * ToRankine(AirTemperature)))
  End Function


  Function ToFahrenheit(ByVal Celsius As Double) As Double
    Return Celsius * 1.8 + 32
  End Function

  Function ToCelsius(ByVal Fahrenheit As Double) As Double
    Return (Fahrenheit - 32) / 1.8
  End Function

  Function ToRadians(ByVal Degrees As Double) As Double
    Return Degrees * π / 180
  End Function

  Function ToDegrees(ByVal Radians As Double) As Double
    Return Radians * 180 / π
  End Function

  Function ToMeters(ByVal Feet As Double) As Double
    Return Feet * 0.304800609601219
  End Function

  Function ToFeet(ByVal Meters As Double) As Double
    Return Meters / 0.304800609601219
  End Function

  Function ToInches(ByVal Millimeters As Double) As Double
    Return Millimeters / 25.4
  End Function

  Function ToMJperM2(ByVal Langleys) As Double
    Return 0.04184 * Langleys
  End Function

  Function ToKelvin(ByVal Celsius) As Double
    Return 273.16 + Celsius
  End Function

  Function ToRankine(ByVal Fahrenheit) As Double
    Return Fahrenheit + 459.67
  End Function

  Function ToLangleysPerHour(ByVal WattsPerSquaredMeter As Double) As Double
    Return WattsPerSquaredMeter / 11.63
  End Function

  Function ToMilesPerHour(ByVal MetersPerSecond As Double) As Double
    Return MetersPerSecond / 0.44704
  End Function

  Function ToMetersPerSecond(ByVal MilesPerHour As Double) As Double
    Return MilesPerHour * 0.44704
  End Function

  ''' <summary>
  ''' Container for geocentric sun angles and distance at a given instant in time.
  ''' </summary>
  Class SolarPosition
    ''' <summary>Geocentric Sun Right Ascension (radians)</summary>
    Public α As Double

    ''' <summary>Geocentric Sun Declination (radians)</summary>
    Public δ As Double

    ''' <summary>Earth-Sun Distance (astronomical units)</summary>
    Public R As Double

    ''' <summary>Apparent Sidereal Time at Greenwich (radians)</summary>
    Public ν As Double

    ''' <summary>
    '''
    ''' </summary>
    ''' <param name="R">Earth-Sun Distance (astronomical units)</param>
    ''' <param name="δ">Geocentric Sun Declination (radians)</param>
    ''' <param name="α">Geocentric Sun Right Ascension (radians)</param>
    ''' <param name="ν">Apparent Sidereal Time at Greenwich (radians)</param>
    ''' <remarks></remarks>
    Sub New(ByVal R, ByVal δ, ByVal α, ByVal ν)
      Me.R = R
      Me.δ = δ
      Me.α = α
      Me.ν = ν
    End Sub
  End Class

  Function Limit(ByRef Value As Double, ByRef Lower As Double, ByRef Upper As Double) As Double
    If Value < Lower Then
      Return Lower
    ElseIf Value > Upper Then
      Return Upper
    Else
      Return Value
    End If
  End Function

  Enum ReferenceET
    LongReference
    ShortReference
  End Enum

  Enum EarthHeliocentricType
    Longitude
    Latitude
    Radius
  End Enum

  Enum NutationType
    Longitude
    Obliquity
  End Enum


  Public Sub Main(args() As string)
    Console.WriteLine("Hello, World!")

    'Dim Output As Double = AdjustAirPressure(85.88, 88.57, -6.24)
    'Console.WriteLine(Output)

    'Dim Output As Double = CalculateWindSpeedAdjustmentFactor(ToFeet(10))
    'Console.WriteLine(Output)

    ' 'date, lon, lat, elev, slope, azimuth, temp, pressure, solar_position
    ' Dim TestDate As DateTime = New DateTime(2003, 10, 17, 19, 30, 30, 0, DateTimeKind.Utc)
    ' 'Dim TestDate As DateTime = New DateTime(2017, 7, 1, 20, 0, 0, 0, DateTimeKind.Utc)
    ' Console.WriteLine(TestDate)
    ' Dim TestSP As SolarPosition = CalculateSolarPosition(TestDate)
    ' Console.WriteLine(TestSP.R)
    ' Console.WriteLine(TestSP.δ)
    ' Console.WriteLine(TestSP.α)
    ' Console.WriteLine(TestSP.ν)

    ' ' date, lon, lat, elev, slope, azimuth, temp, pressure, solar_position
    ' Dim TestDate As DateTime = New DateTime(2017, 7, 1, 18, 0, 0, 0, DateTimeKind.Utc)
    ' 'Dim TestDate As DateTime = New DateTime(2017, 7, 1, 20, 0, 0, 0, DateTimeKind.Utc)
    ' Dim Test_SP As SolarPosition = CalculateSolarPosition(TestDate)
    ' Dim Output As Double = CalculateInstantaneousRa(
    '   -112.688, 39.313, 4584, 0.04, 275.6, 86.27, 85.979, Test_SP
    ' )
    ' Console.WriteLine(Output)

    ' Dim TimeStep = 15 'minutes
    ' Dim HourlySteps As Integer = 60 / TimeStep
    ' Dim SolarPositions(1440 / TimeStep) As SolarPosition
    ' Dim RecordDate As DateTime = New DateTime(2017, 7, 1, 13, 0, 0, 0, DateTimeKind.Utc)
    ' For I as Integer = 0 To SolarPositions.Length - 1
    '     SolarPositions(I) = CalculateSolarPosition(RecordDate.AddMinutes(I * TimeStep - 30))
    ' Next

    ' ' lon, lat, elev, slope, azimuth, temp, pressure, solar_position
    ' Dim TestDate As DateTime = New DateTime(2017, 7, 1, 0, 0, 0, 0, DateTimeKind.Utc)
    ' Dim Test_SP As SolarPosition = CalculateSolarPosition(TestDate)
    ' ' Test point 1
    ' Dim Output1 As Double = CalculateInstantaneousRa(
    '   -113.42677351013316, 41.81679884002846, 5669.700337414034, 1.691116213798523, 189.4030303955078,
    '   72.10305476625949, 83.11632613999164, Test_SP
    ' )
    ' Console.WriteLine(Output1)
    ' ' Test point 2
    ' Dim Output2 As Double = CalculateInstantaneousRa(
    '   -113.41373207120276, 41.812179783631834, 5593.457918405696, 2.260589599609375, 164.3304443359375,
    '   72.46959210387394, 83.32195663998966, Test_SP
    ' )
    ' Console.WriteLine(Output2)

    ' ' RsSource, RaSource, RaDestination, Slope
    ' Dim Output As Double = TranslateRs(0.42061625, 0.0567711, 0.38273549857, 2.15562152863)
    ' Console.WriteLine(Output)
    ' Dim Output As Double = TranslateRs(0.42061625, 0.0567711, 0, 2.26058959961)
    ' Console.WriteLine(Output)

    ' ' elev, temp, rh, wind, rs, ra, zw, type, pres
    ' Dim Output As Double = CalculateHourlyASCEReferenceET(
    '   4584, 86.27, 11.88, 0.786, 75, 105, 2 / 0.3048, ReferenceET.ShortReference, 85.979
    ' )
    ' Console.WriteLine(Output)

    Dim MinDate As DateTime = New DateTime(2017, 7, 1, 13, 0, 0, 0, DateTimeKind.Utc)
    Dim MaxDate As DateTime = New DateTime(2017, 7, 2, 12, 0, 0, 0, DateTimeKind.Utc)
    Dim Output As Double = CalculateReferenceEvapotranspirationNLDAS_2A(
      MinDate, MaxDate,
      -112.6854126892207, 39.31432060447445, 4584.507531242295, 0.07245257496833801, 239.5464324951172,
      69.08371356284556, 86.10269479690949,
      -112.6875, 39.3125, 4613.520889383954, 0.038626529276371, 263.5475158691406,
      69.098, 86.10418, 9.678417884780739, 15.147373521610197, -3.154080171796707, -1.7671796707229779
    )
    Console.WriteLine(Output)

  End Sub

End Module