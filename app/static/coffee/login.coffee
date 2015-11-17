class Login
  view: -> m.request(
    {
        method : 'GET'
        url : "http://127.0.0.1:5000"
        headers: {'Access-Control-Allow-Origin: *'},
    })

login = new Login()
