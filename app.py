from flask import Flask, request, jsonify
from flask_cors import CORS
import sympy as sp
import re
import sys
import threading
import time
import requests

app = Flask(__name__)
CORS(app)

latest_response = None

class ODE:
  def __init__(self, odeEquation, conditionList):
    self.a = sp.symbols("a")
    self.b = sp.symbols("b")
    self.c = sp.symbols("c")
    self.d = sp.symbols("d")
    self.e = sp.symbols("e")
    self.f = sp.symbols("f")
    self.g = sp.symbols("g")
    self.h = sp.symbols("h")
    self.i = sp.symbols("i")
    self.j = sp.symbols("j")
    self.k = sp.symbols("k")
    self.l = sp.symbols("l")
    self.m = sp.symbols("m")
    self.n = sp.symbols("n")
    self.o = sp.symbols("o")
    self.p = sp.symbols("p")
    self.q = sp.symbols("q")
    self.r = sp.symbols("r")
    self.s = sp.symbols("s")
    self.t = sp.symbols("t")
    self.u = sp.symbols("u")
    self.v = sp.symbols("v")
    self.w = sp.symbols("w")
    self.x = sp.symbols("x")
    self.z = sp.symbols("z")
    self.y = sp.Function("y")(self.x)
    self.order = 0
    self.isHomogenous = True
    self.odeEquation = odeEquation
    self.odeEquation = self.getEquation()
    self.conditionList = conditionList
    self.errorString = ""
    self.isEquation = 1
    self.odeSolution = 0
    self.conditions = {}
    if type(self.odeEquation) == str:
      self.errorString = self.odeEquation
      self.isEquation = 0

  def getConditionList(self):
    return self.conditionList

  def isEquationFormed(self):
    return self.isEquation

  def getError(self):
    return self.errorString

  def isTermEnd(self, odeEquationTerm, i):
    if odeEquationTerm[i - 1] == "D" or odeEquationTerm[i - 1] == "y":
      return True
    for j in range(2, i):
      if odeEquationTerm[i -
                         j] == "^" and (odeEquationTerm[i - j - 1] == "D"
                                        or odeEquationTerm[i - j - 1] == "y"):
        return True
      if odeEquationTerm[i - j] == "+" or odeEquationTerm[i - j] == "-":
        return False
    return False

  def extractYfromD(self, dCoEfficient):
    self.isHomogenous = False
    if "*y^" in dCoEfficient:
      parts = dCoEfficient.split("*y^")
      YCoEfficient = parts[0]
      power = parts[1]
      return YCoEfficient, power
    elif "y^" in dCoEfficient:
      parts = dCoEfficient.split("y^")
      YCoEfficient = parts[0]
      if YCoEfficient == "-" or YCoEfficient == "+" or YCoEfficient == "":
        YCoEfficient += "1"
      power = parts[1]
      return YCoEfficient, power
    elif "*y" in dCoEfficient:
      parts = dCoEfficient.split("*y")
      YCoEfficient = parts[0]
      power = 1
      return YCoEfficient, power
    elif "y" in dCoEfficient:
      parts = dCoEfficient.split("y")
      YCoEfficient = parts[0]
      if YCoEfficient == "-" or YCoEfficient == "+" or YCoEfficient == "":
        YCoEfficient += "1"
      power = 1
      return YCoEfficient, power

  def extractTerms(self, odeEquationTerm):
    if odeEquationTerm == "0" or odeEquationTerm == "+0" or odeEquationTerm == "-0":
      return [0], [0], [], ["0"], ["0"]
    powerListOfD = []
    powerListOfYInD = []
    powerListOfY = []
    coEfficientListOfD = []
    coEfficientListOfY = []
    term = ""
    YCoEfficient = ""
    dCoEfficient = ""
    yPower = 0
    isBracketsFullClosed = 0
    for i in range(0, len(odeEquationTerm)):
      char = odeEquationTerm[i]
      if char == "(":
        isBracketsFullClosed += 1
      if char == ")":
        isBracketsFullClosed -= 1
      if (i and char in ("+", "-") and odeEquationTerm[i - 1] != "^" and
          (self.isTermEnd(odeEquationTerm, i) or isBracketsFullClosed == 0)):
        if "D" not in term and "y" not in term:
          dCoEfficient = term
          power = -sys.maxsize
          yPower = 0
          powerListOfYInD.append(float(yPower))
          powerListOfD.append(int(power))
        elif term:
          if "*D^" in term:
            parts = term.split("*D^")
            dCoEfficient = parts[0]
            if "y" in dCoEfficient:
              dCoEfficient, yPower = self.extractYfromD(dCoEfficient)
            else:
              yPower = 0
            powerListOfYInD.append(float(yPower))
            power = parts[1]
            powerListOfD.append(int(power))
          elif "D^" in term:
            parts = term.split("D^")
            dCoEfficient = parts[0]
            if (dCoEfficient == "-" or dCoEfficient == "+"
                or dCoEfficient == ""):
              dCoEfficient += "1"
            if "y" in dCoEfficient:
              dCoEfficient, yPower = self.extractYfromD(dCoEfficient)
            else:
              yPower = 0
            powerListOfYInD.append(float(yPower))
            power = parts[1]
            powerListOfD.append(int(power))
          elif "*D" in term:
            parts = term.split("*D")
            dCoEfficient = parts[0]
            if "y" in dCoEfficient:
              dCoEfficient, yPower = self.extractYfromD(dCoEfficient)
            else:
              yPower = 0
            powerListOfYInD.append(float(yPower))
            power = 1
            powerListOfD.append(int(power))
          elif "D" in term:
            parts = term.split("D")
            dCoEfficient = parts[0]
            if (dCoEfficient == "-" or dCoEfficient == "+"
                or dCoEfficient == ""):
              dCoEfficient += "1"
            if "y" in dCoEfficient:
              dCoEfficient, yPower = self.extractYfromD(dCoEfficient)
            else:
              yPower = 0
            powerListOfYInD.append(float(yPower))
            power = 1
            powerListOfD.append(int(power))
          elif "*y^" in term:
            parts = term.split("*y^")
            YCoEfficient = parts[0]
            power = parts[1]
            powerListOfY.append(float(power))
          elif "y^" in term:
            parts = term.split("y^")
            YCoEfficient = parts[0]
            if (YCoEfficient == "-" or YCoEfficient == "+"
                or YCoEfficient == ""):
              YCoEfficient += "1"
            power = parts[1]
            powerListOfY.append(float(power))
          elif "*y" in term:
            parts = term.split("*y")
            YCoEfficient = parts[0]
            power = 1
            powerListOfY.append(power)
          elif "y" in term:
            parts = term.split("y")
            YCoEfficient = parts[0]
            if (YCoEfficient == "-" or YCoEfficient == "+"
                or YCoEfficient == ""):
              YCoEfficient += "1"
            power = 1
            powerListOfY.append(power)
        if dCoEfficient:
          coEfficientListOfD.append(dCoEfficient)
        if YCoEfficient:
          coEfficientListOfY.append(YCoEfficient)
        term = char
        YCoEfficient = ""
        dCoEfficient = ""
      else:
        term += char
    if "D" not in term and "y" not in term:
      dCoEfficient = term
      power = -sys.maxsize
      yPower = 0
      powerListOfYInD.append(float(yPower))
      powerListOfD.append(int(power))
    elif term:
      if "*D^" in term:
        parts = term.split("*D^")
        dCoEfficient = parts[0]
        if "y" in dCoEfficient:
          dCoEfficient, yPower = self.extractYfromD(dCoEfficient)
        else:
          yPower = 0
        powerListOfYInD.append(float(yPower))
        power = parts[1]
        powerListOfD.append(int(power))
      elif "D^" in term:
        parts = term.split("D^")
        dCoEfficient = parts[0]
        if dCoEfficient == "-" or dCoEfficient == "+" or dCoEfficient == "":
          dCoEfficient += "1"
        if "y" in dCoEfficient:
          dCoEfficient, yPower = self.extractYfromD(dCoEfficient)
        else:
          yPower = 0
        powerListOfYInD.append(float(yPower))
        power = parts[1]
        powerListOfD.append(int(power))
      elif "*D" in term:
        parts = term.split("*D")
        dCoEfficient = parts[0]
        if "y" in dCoEfficient:
          dCoEfficient, yPower = self.extractYfromD(dCoEfficient)
        else:
          yPower = 0
        powerListOfYInD.append(float(yPower))
        power = 1
        powerListOfD.append(int(power))
      elif "D" in term:
        parts = term.split("D")
        dCoEfficient = parts[0]
        if dCoEfficient == "-" or dCoEfficient == "+" or dCoEfficient == "":
          dCoEfficient += "1"
        if "y" in dCoEfficient:
          dCoEfficient, yPower = self.extractYfromD(dCoEfficient)
        else:
          yPower = 0
        powerListOfYInD.append(float(yPower))
        power = 1
        powerListOfD.append(int(power))
      elif "*y^" in term:
        parts = term.split("*y^")
        YCoEfficient = parts[0]
        power = parts[1]
        powerListOfY.append(float(power))
      elif "y^" in term:
        parts = term.split("y^")
        YCoEfficient = parts[0]
        if YCoEfficient == "-" or YCoEfficient == "+" or YCoEfficient == "":
          YCoEfficient += "1"
        power = parts[1]
        powerListOfY.append(float(power))
      elif "*y" in term:
        parts = term.split("*y")
        YCoEfficient = parts[0]
        power = 1
        powerListOfY.append(power)
      elif "y" in term:
        parts = term.split("y")
        YCoEfficient = parts[0]
        if YCoEfficient == "-" or YCoEfficient == "+" or YCoEfficient == "":
          YCoEfficient += "1"
        power = 1
        powerListOfY.append(power)
    if dCoEfficient:
      coEfficientListOfD.append(dCoEfficient)
    if YCoEfficient:
      coEfficientListOfY.append(YCoEfficient)
    return (
        powerListOfD,
        powerListOfYInD,
        powerListOfY,
        coEfficientListOfD,
        coEfficientListOfY,
    )

  def convertEquation(self, odeEquationTerm):
    odeEquationTerm = odeEquationTerm.replace("e^", "exp")
    odeEquationTerm = odeEquationTerm.replace("root", "sqrt")
    odeEquationTerm = odeEquationTerm.replace("d", "D")
    odeEquationTerm = odeEquationTerm.replace("X", "x")
    odeEquationTerm = odeEquationTerm.replace("Y", "y")
    odeEquationTerm = odeEquationTerm.replace("(D)", "D")
    odeEquationTerm = odeEquationTerm.replace("(y)", "y")
    odeEquationTerm = odeEquationTerm.replace("D*", "D")
    odeEquationTerm = odeEquationTerm.replace("y*", "y")
    pattern = r"(D\^-?\d*)y"
    odeEquationTerm = re.sub(pattern, r"\1", odeEquationTerm)
    pattern = r"(D\^\+?\d*)y"
    odeEquationTerm = re.sub(pattern, r"\1", odeEquationTerm)
    return odeEquationTerm

  def getEquation(self):
    errorMessages = {
        "y/": "Don't divide the term 'y', simplify the equation even more.",
        "y)": "Don't wrap 'y' with ().",
        "D)": "Don't wrap 'D' with ().",
        "(D": "Don't wrap 'D^n' with ().",
        "D^(": "Don't wrap D's power with ().",
        "y^(": "Don't wrap y's power with ().",
    }
    try:
      self.odeEquation = self.convertEquation(self.odeEquation.replace(
          " ", ""))
      for condition, errorMessage in errorMessages.items():
        if condition in self.odeEquation:
          return errorMessage
      LhsOde, RhsOde = self.odeEquation.split("=")
    except:
      LhsOde = self.odeEquation
      RhsOde = "0"
    (
        LhsPowerListOfD,
        LhsPowerListOfyInD,
        LhsPowerListOfY,
        LhsCoEfficientListOfD,
        LhsCoEfficientListOfY,
    ) = self.extractTerms(LhsOde)
    (
        RhsPowerListOfD,
        RhsPowerListOfyInD,
        RhsPowerListOfY,
        RhsCoEfficientListOfD,
        RhsCoEfficientListOfY,
    ) = self.extractTerms(RhsOde)
    if not LhsPowerListOfD and not RhsPowerListOfD:
      return "Order has to be greater than 0."
    if len(LhsPowerListOfY) and max(LhsPowerListOfY) != 1:
      self.isHomogenous = False
    if len(RhsPowerListOfY) and max(RhsPowerListOfY) != 1:
      self.isHomogenous = False
    for num in LhsPowerListOfD:
      if num < 0 and num != -sys.maxsize:
        return "Order has to be greater than 0."
      if num == -sys.maxsize:
        self.isHomogenous = False
    for num in RhsPowerListOfD:
      if num < 0 and num != -sys.maxsize:
        return "Order has to be greater than 0."
      if num == -sys.maxsize:
        self.isHomogenous = False
    if "x" in self.odeEquation:
      self.isHomogenous = False
    if (len(LhsPowerListOfD) == 1 and RhsOde == "0"
        and len(LhsCoEfficientListOfY) == 0):
      self.isHomogenous = True
    if len(LhsPowerListOfD) and len(RhsPowerListOfD):
      self.order = max(max(LhsPowerListOfD), max(RhsPowerListOfD))
    elif len(LhsPowerListOfD):
      self.order = max(LhsPowerListOfD)
    else:
      self.order = max(RhsPowerListOfD)
    if self.order == 0:
      return "Order has to be greater than 0."
    LhsEquation = 0
    for i, item in enumerate(LhsPowerListOfD):
      if item == -sys.maxsize:
        LhsEquation += sp.sympify(LhsCoEfficientListOfD[i])
      elif isinstance(LhsPowerListOfyInD[i],
                      float) and LhsPowerListOfyInD[i] == int(
                          LhsPowerListOfyInD[i]):
        LhsEquation += (self.y.diff(self.x, item) *
                        (self.y**int(LhsPowerListOfyInD[i])) *
                        sp.sympify(LhsCoEfficientListOfD[i]))
      elif isinstance(
          LhsPowerListOfyInD[i],
          float) and LhsPowerListOfyInD[i] != int(LhsPowerListOfyInD[i]):
        LhsEquation += (self.y.diff(self.x, item) *
                        (self.y**LhsPowerListOfyInD[i]) *
                        sp.sympify(LhsCoEfficientListOfD[i]))
    for i, item in enumerate(LhsPowerListOfY):
      if isinstance(item, float) and item != int(item):
        LhsEquation += (self.y**item) * sp.sympify(LhsCoEfficientListOfY[i])
      elif isinstance(item, float) and item == int(item):
        LhsEquation += (self.y**int(item)) * sp.sympify(
            LhsCoEfficientListOfY[i])
      else:
        LhsEquation += self.y * sp.sympify(LhsCoEfficientListOfY[i])
    RhsEquation = 0
    for i, item in enumerate(RhsPowerListOfD):
      if item == -sys.maxsize:
        RhsEquation += sp.sympify(RhsCoEfficientListOfD[i])
      elif isinstance(RhsPowerListOfyInD[i],
                      float) and RhsPowerListOfyInD[i] == int(
                          RhsPowerListOfyInD[i]):
        RhsEquation += (self.y.diff(self.x, item) *
                        (self.y**int(RhsPowerListOfyInD[i])) *
                        sp.sympify(RhsCoEfficientListOfD[i]))
      elif isinstance(
          RhsPowerListOfyInD[i],
          float) and RhsPowerListOfyInD[i] != int(RhsPowerListOfyInD[i]):
        RhsEquation += (self.y.diff(self.x, item) *
                        (self.y**RhsPowerListOfyInD[i]) *
                        sp.sympify(RhsCoEfficientListOfD[i]))
    for i, item in enumerate(RhsPowerListOfY):
      if isinstance(item, float) and item != int(item):
        RhsEquation += (self.y**item) * sp.sympify(RhsCoEfficientListOfY[i])
      elif isinstance(item, float) and item == int(item):
        RhsEquation += (self.y**int(item)) * sp.sympify(
            RhsCoEfficientListOfY[i])
      else:
        RhsEquation += self.y * sp.sympify(RhsCoEfficientListOfY[i])
    equation = sp.Eq(LhsEquation, RhsEquation)
    return equation

  def getHomogeneity(self):
    return self.isHomogenous

  def getOrder(self):
    return self.order

  def getEnteredOde(self):
    return self.odeEquation

  def isGeneralSolution(self):
    try:
      self.odeSolution = sp.dsolve(self.odeEquation)
      return 1
    except:
      return 0

  def getGeneralSolution(self):
    return self.odeSolution

  def getOdeConditions(self):
    return self.conditions

  def getParticularSolution(self):
    for i in range(0, int(len(self.conditionList))):
      for j in range(0, int(len(self.conditionList[i]) / 2)):
        x = self.conditionList[i][2 * j]
        y = self.conditionList[i][2 * j + 1]
        if i:
          self.conditions[self.y.diff(self.x, i).subs(self.x, x)] = y
        else:
          self.conditions[self.y.subs(self.x, x)] = y
    try:
      particularSolution = sp.dsolve(self.odeEquation, ics=self.conditions)
    except:
      return 0
    return particularSolution


@app.route("/odeSolve", methods=["POST"])
def mainFun():
  try:
    data = request.json
    dataEquation = data.get("equation", "")
    conditionList = data.get("inCon", [])

    print("\n", conditionList, "\n")
    print("\n" + dataEquation + "\n")

    odeSolver = ODE(dataEquation, conditionList)
    enteredEq = 0
    generalEq = 0
    particularEq = 0

    flag = "Y" if len(odeSolver.getConditionList()) != 0 else "n"

    if flag == "Y":
      particularSolution = odeSolver.getParticularSolution()
      try:
        if particularSolution:
          print("\nCan't Solve, enter conditions carefully.\n")
          return ["Can't Solve, enter conditions carefully."], 200
      except:
        particularEq = str(sp.latex(particularSolution))
        print("\n" + sp.pretty(particularSolution, use_unicode=True) + "\n")
        return [particularEq], 200

    if odeSolver.isEquationFormed():
      enteredEq = odeSolver.getEnteredOde()
      print("\n" + sp.pretty(enteredEq, use_unicode=True) + "\n")
      enteredEq = str(sp.latex(enteredEq))
    else:
      string = odeSolver.getError()
      return [string], 200

    if odeSolver.isGeneralSolution():
      generalSolution = odeSolver.getGeneralSolution()
      print("\n" + sp.pretty(generalSolution, use_unicode=True) + "\n")
      generalEq = str(sp.latex(generalSolution))
    else:
      print("\nCan't Solve.\n")
      return ["Can't Solve."], 200

    if odeSolver.getHomogeneity():
      print(f"\nIt's a homogenous equation. order {odeSolver.getOrder()}\n")
      return [
          enteredEq,
          generalEq,
          f"\nIt's a homogenous equation. order {odeSolver.getOrder()}",
          dataEquation,
      ]
    else:
      print(
          f"\nIt's not a homogenous equation. order {odeSolver.getOrder()}\n")
      return [
          enteredEq,
          generalEq,
          f"\nIt's not a homogenous equation. order {odeSolver.getOrder()}",
          dataEquation,
      ]
  except:
    print("\nAn Unexpected error occured.\n")
    return ["An Unexpected error occured."], 200


def fetch_data():
    global latest_response
    while True:
        try:
            response = requests.get("https://omrevalserver.onrender.com/")
            latest_response = response
            print(latest_response)
        except Exception as e:
            print(f"Error fetching data: {e}")
        time.sleep(1)


@app.route("/")
def get_data():
    if latest_response is not None:
        return jsonify({"message": "omrElite"}), 200
    else:
        return jsonify({"message": "No data available"}), 404
