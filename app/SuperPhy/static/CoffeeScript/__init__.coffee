###
    Events
###

#Routes are a set of endpoints that the user can navigate to without losing
#the SPA context. They are added at the Page level.

Routes = {}
Routes.add = (route, class_, args) ->
    args = args || {}
    Routes.routes ?= {}
    if route.substring(0,1) is '/'
        route = route.substring(1)
    Routes.routes["/#{route}"] = view: -> m.component(class_, args)
