import datetime as dt

import pandas as pd
import requests
from flask import Flask, request
from flask_cors import CORS
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.linear_model import LinearRegression
from sklearn.metrics import mean_squared_error
from sklearn.model_selection import train_test_split

app = Flask(__name__)
CORS(app)

sun_shine_hour = [
    [7.5, 8, 8.25, 8, 7.5, 4.7, 4.4, 4.5, 5.2, 7, 7.6, 7.9],
    [6.6, 7.8, 7.8, 8, 7.2, 5.4, 4.2, 5.2, 5.7, 6.4, 7.3, 7.4],
    [8.2, 8.6, 8.7, 8.6, 7.6, 4.9, 4.5, 4.9, 5.8, 7.4, 8, 8.1],
    [6.7, 7.5, 7.7, 7.9, 6.7, 5.8, 5.2, 5.5, 5.5, 6.6, 7.4, 7.3],
    [6.8, 7.5, 8, 8.2, 7.8, 4.8, 4.7, 5.1, 4.8, 6.5, 7.5, 6.7],
    [6.1, 7.5, 8, 7, 6.8, 4.8, 5, 5.3, 5.3, 6.7, 7.6, 6],
    [7.1, 7.9, 8, 7.8, 7.3, 5.6, 4.5, 5.4, 5.5, 7.4, 7.7, 7.1],
    [6.3, 7.7, 7.4, 7.5, 7.3, 4.7, 3.9, 4.2, 4.4, 6.8, 7.1, 6.3],
    [6.1, 7.2, 7.9, 6.9, 6.7, 4.2, 4.4, 4.1, 4.2, 7.2, 7.7, 6.7],
    [6.8, 8.2, 8.3, 8.2, 7.9, 6.8, 4.4, 5.2, 6.4, 7.4, 7.8, 7.7],
    [7.8, 8.2, 8.1, 8, 6.2, 4.5, 3.8, 5.7, 5.9, 6.3, 7.9, 8],
    [7.7, 7.8, 8, 7.4, 6.5, 4.5, 3.8, 4.8, 5, 7.1, 8.1, 8],
    [6.4, 7.8, 7.9, 7.8, 7.1, 5.5, 4.8, 5.6, 5, 7.1, 8.2, 6.8],
    [7, 8.1, 8.3, 8.8, 8.1, 5.5, 4.3, 4.7, 5, 6.8, 7.7, 7.9],
]


def bengali_converter(pre):
    pre_int = int(pre) - 1
    pre_float = pre - pre_int
    bangla_list = [
        ["Poush", "Magh"],
        ["Magh", "Falgun"],
        ["Falgun", "chaitra"],
        ["chaitro", "Boisakh"],
        ["Boisakh", "Joistho"],
        ["Joistho", "Ashar"],
        ["Ashar", "shraban"],
        ["shraban", "bhadro"],
        ["bhadro", "ashwin"],
        ["ashwin", "kartik"],
        ["kartik", "agrahayan"],
        ["agrahayan", "Poush"],
    ]
    bangla_month_list = bangla_list[pre_int]
    if pre_float > 0.5:
        bangla_month = bangla_month_list[1]
        if pre_float > 0.75:
            bangla_week = 2
        else:
            bangla_week = 1
    else:
        bangla_month = bangla_month_list[0]
        if pre_float > 0.25:
            bangla_week = 2
        else:
            bangla_week = 1
    return bangla_month, bangla_week


def get_month_and_week(number):
    number_int = int(number) - 1
    number_float = number - number_int
    if number_float > 0.5:
        if number_float > 0.75:
            week = 4
        else:
            week = 3
    else:
        if number_float > 0.25:
            week = 2
        else:
            week = 1
    return number_int, week


def model2(inputs, outputs, input_list):
    INFO = "Agriqo(slider2).csv"
    data = pd.read_csv(INFO)

    data_without_duplicates = data.drop_duplicates()

    model = LinearRegression()
    model.fit(data_without_duplicates[inputs], data_without_duplicates[outputs])
    predicted_values = model.predict(data_without_duplicates[inputs])
    mse = mean_squared_error(data_without_duplicates[outputs], predicted_values)
    pre = model.predict([input_list])
    pre = pre[0][0]
    bangla_month, bangla_week = bengali_converter(pre)
    english_month, eng_week = get_month_and_week(pre)

    return pre, bangla_month, bangla_week, english_month, eng_week


def weather_api(agricultural_zone, city):
    api_key = "544ee17b2e3eafb5c5a4ec602c2a19b4"
    base_url = "https://api.openweathermap.org/data/2.5/forecast?q="

    Cc = "BD"

    url = base_url + city + "," + Cc + "&appid=" + api_key

    data = requests.get(url).json()

    temperatures = [item["main"]["temp"] for item in data["list"]]
    average_temperature = sum(temperatures) / len(temperatures)

    humidities = [item["main"]["humidity"] for item in data["list"]]
    humidities.sort(reverse=True)
    num_to_avg = min(15, len(humidities))
    average_humidity = sum(humidities[:num_to_avg]) / (num_to_avg - 4.9)

    wind_speeds = [item["wind"]["speed"] for item in data["list"]]
    average_wind_speed = sum(wind_speeds) / len(wind_speeds)

    wind_directions = [
        item["wind"]["deg"] for item in data["list"] if "deg" in item["wind"]
    ]
    average_wind_direction = sum(wind_directions) / len(wind_directions)

    rains = [item["rain"]["3h"] for item in data["list"] if "rain" in item]
    if rains:
        average_rainfall = sum(rains)
    else:
        average_rainfall = 0

    current_date = dt.datetime.today().day
    current_monthe = dt.datetime.today().month
    current_month = current_monthe

    if current_date >= 1 and current_date <= 7:
        current_month += 0.0
    elif current_date >= 8 and current_date <= 14:
        current_month += 0.25
    elif current_date >= 15 and current_date <= 21:
        current_month += 0.5
    elif current_date >= 8 and current_date <= 7:
        current_month += 0.75

    sun_light_hour = sun_shine_hour[agricultural_zone][current_monthe] * 7

    return (
        current_month,
        average_rainfall,
        average_temperature,
        average_humidity,
        sun_light_hour,
        average_wind_direction,
        average_wind_speed,
    )


def model1(inputs, outputs, input_list):
    INFO = "Agriqo(slider1).csv"
    data = pd.read_csv(INFO)
    X_train, X_test, Y_train, Y_test = train_test_split(
        data[inputs], data[outputs], test_size=0.2, random_state=42
    )
    model = HistGradientBoostingClassifier(random_state=42)
    model.fit(X_train, Y_train)
    probabilites = model.predict_proba([input_list])[0]
    top_3 = probabilites.argsort()[-3:][::-1]
    top_3_classes = model.classes_[top_3]
    top_3_probabilities = probabilites[top_3]

    return top_3_classes, top_3_probabilities


def ai(inputs, manual_input):
    if not manual_input:
        (
            current_month,
            average_rainfall,
            average_temperature,
            average_humidity,
            sun_light_hour,
            average_wind_direction,
            average_wind_speed,
        ) = weather_api(inputs[0], inputs[1])
        input_list = [
            inputs[0],
            current_month,
            average_rainfall,
            average_temperature,
            average_humidity,
            sun_light_hour,
            average_wind_direction,
            average_wind_speed,
        ]
        output_class, output = model1(
            [
                "Agricultural zone",
                "month(chara)",
                "rainfall (mm)",
                "temperature(avg)",
                "humidity(avg)",
                "sunlight(hour)",
                "direction of wind(deg)",
                "velocity of wind(km/h)",
            ],
            ["label"],
            input_list,
        )
    elif manual_input:
        current_date = dt.datetime.today().day
        current_monthe = dt.datetime.today().month
        current_month = current_monthe

        if current_date >= 1 and current_date <= 7:
            current_month += 0.0
        elif current_date >= 8 and current_date <= 14:
            current_month += 0.25
        elif current_date >= 15 and current_date <= 21:
            current_month += 0.5
        elif current_date >= 8 and current_date <= 7:
            current_month += 0.75
        input_list = [
            inputs[0],
            current_month,
            inputs[2],
            inputs[3],
            inputs[4],
            inputs[5],
            inputs[6],
            inputs[7],
        ]
        output_class, output = model1(
            [
                "Agricultural zone",
                "month(chara)",
                "rainfall (mm)",
                "temperature(avg)",
                "humidity(avg)",
                "sunlight(hour)",
                "direction of wind(deg)",
                "velocity of wind(km/h)",
            ],
            ["label"],
            input_list,
        )

    output_str = "["

    for i, j in zip(output, output_class):
        if i > 0.3:
            output_str += '"' + j + '",'

    output_list = list(output_str)
    output_list[len(output_str) - 1] = "]"
    final_output = "".join(output_list)

    return final_output


@app.route("/cropRecomAI", methods=["POST"])
def handle_data_Ai():
    try:
        data = request.json

        inputs = [data["zone"], data["district"]]
        manual_input = False

        return ai(inputs, manual_input)
    except:
        return "An Unexpected error occured.", 200


@app.route("/cropRecomCustom", methods=["POST"])
def handle_data_custom():
    try:
        data = request.json

        inputs = [
            data["zone"],
            data["district"],
            float(data["rainfall"]),
            float(data["temperature"]),
            float(data["humidity"]),
            float(data["sunshine"]),
            float(data["windDirection"]),
            float(data["windVelocity"]),
        ]
        manual_input = True

        return ai(inputs, manual_input)
    except:
        return "An Unexpected error occured.", 200


@app.route("/timeRecomAI", methods=["POST"])
def handle_data_time():
    try:
        data = request.json

        input_list = [data["zone"], data["crop"]]
        months = [
            "January",
            "February",
            "March",
            "April",
            "May",
            "June",
            "July",
            "August",
            "September",
            "October",
            "November",
            "December",
        ]
        pre, bangla_month, bangla_week, english_month, eng_week = model2(
            ["Agricultural zone", "label count"], ["month(chara)"], input_list
        )

        string = f"week {eng_week} of {months[english_month]}<br/>week {bangla_week} of {bangla_month}"

        return string
    except:
        return "An Unexpected error occured.", 200


@app.route("/")
def handle_data():
    return "Hello"
