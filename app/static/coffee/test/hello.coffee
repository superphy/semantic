
class Hello
    urlhello= "http://10.139.14.121:5000/mithril/query"

    controller =
        getData: -> m.request({
          method: "GET",
          url: urlhello
          background: true,
          initialValue: []
        })

    model=
      flowers: ["roses","violets"]
      data: controller.getData()

    view: () ->
        [
            home.view()
            m("div", ["Getting the first 3 triples in the triplestore."])
            m("div", model.flowers)
            m("div", JSON.stringify(model.data))
            #m("div",[(data.results)])
#            m("p.home-beta", {class:'text-center'}, [
#                JSON.stringify(users)
#            ])
        ]

