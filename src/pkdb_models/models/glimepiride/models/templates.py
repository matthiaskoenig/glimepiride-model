"""Template definition."""
from sbmlutils.factory import *
from datetime import datetime

creators = [
    Creator(
        familyName="Elias",
        givenName="Michelle",
        email="michelleelias66@gmail.com",
        organization="Humboldt-University Berlin, Institute for Biology",
    ),
    Creator(
       familyName="König",
       givenName="Matthias",
       email="koenigmx@hu-berlin.de",
       organization="Humboldt-University Berlin, Institute for Theoretical Biology",
       site="https://livermetabolism.com",
    ),
]

terms_of_use = """
    The content of this model has been carefully created in a manual research effort.

    ## Terms of use
    Copyright © {year} <a href="{site}" title="{given_name} {family_name}" target="_blank">{given_name} {family_name}</a>.

    <a rel="license" href="http://creativecommons.org/licenses/by/4.0/">
    <img alt="Creative Commons License" style="border-width:0" src="https://i.creativecommons.org/l/by/4.0/88x31.png" /></a><br />This work is licensed under a <a rel="license" href="http://creativecommons.org/licenses/by/4.0/">Creative Commons Attribution 4.0 International License</a>.

    Redistribution and use of any part of this model, with or without modification, are permitted provided
    that the following conditions are met:

    - Redistributions of this SBML file must retain the above copyright notice, this list of conditions and the following disclaimer.
    - Redistributions in a different form must reproduce the above copyright notice, this list of conditions
      and the following disclaimer in the documentation and/or other materials provided
      with the distribution.

    This model is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the
    implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
""".format(
    year=datetime.now().year,
    given_name=creators[-1].givenName,
    family_name=creators[-1].familyName,
    site=creators[-1].site,
)


class U(Units):
    """UnitDefinitions."""

    mmole = UnitDefinition("mmole")
    min = UnitDefinition("min")
    mg = UnitDefinition("mg")
    m2 = UnitDefinition("m2", "meter^2")
    mM = UnitDefinition("mM", "mmole/liter")
    mmole_per_min = UnitDefinition("mmole_per_min", "mmole/min")
    mmole_per_min_l = UnitDefinition("mmole_per_min_l", "mmole/min/l")
    l_per_min = UnitDefinition("l_per_min", "l/min")
    per_min = UnitDefinition("per_min", "1/min")
    per_mmole = UnitDefinition("per_mmole", "1/mmole")
    mg_per_s_l = UnitDefinition("mg_per_s_l", "mg/s/l")
    g_per_mole = UnitDefinition("g_per_mole", "g/mole")
    min_per_hr = UnitDefinition("min_per_hr", "min/hour")


model_units = ModelUnits(
    time=U.min,
    extent=U.mmole,
    substance=U.mmole,
    length=U.meter,
    area=U.m2,
    volume=U.liter,
)


