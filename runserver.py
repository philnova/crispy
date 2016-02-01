from crispy import app
app.secret_key = 'super_secret_key'
app.run(host='0.0.0.0', port=8910, debug=True)