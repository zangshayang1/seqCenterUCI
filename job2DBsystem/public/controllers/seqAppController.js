var seqApp = angular.module("seqApp", ["ngRoute"]);

seqApp.controller("seqAppController", ["$scope",
                                       "$http",
                                        function($scope, $http) {
    console.log("job submitted.");
    
    // for use of jobSubmit.html 
    $scope.addSeqJob = function() {
        // to do: set a check point so that only completed forms will go to server
        $http.post("/job2DBsystem", $scope.seq).success(function(response) {
            console.log("new job sent!");    
        });
    };

}]);

seqApp.controller("seqAppAdminController", ["$scope",
                                            "$http",
                                            function($scope, $http) {

    // for use of jobView.html
    // auto load the following when loading the page
    $http.get("/job2DBsystem").success(function(response) {
        $scope.seqJobs = response;
        // $scope.jobSubTime = $scope.seqJobs._id
        console.log($scope.seqJobs); 
        console.log("heard back from server!");
    });    
}]);
