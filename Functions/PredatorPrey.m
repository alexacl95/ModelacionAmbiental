function Y = PredatorPrey(theta, t)
    IC = theta(1:2);
    [~, Y] = ode45(@DifEq,t,IC);

    function dY = DifEq(~,State)
        dY = [theta(2)*State(1)-theta(3)*State(1)*State(2);
              theta(4)*State(1)*State(2) - theta(5)*State(2)];
    end

end