using DataAccess;
using Infrastructure.Models;
using Microsoft.AspNetCore.Mvc;
using Microsoft.AspNetCore.Mvc.RazorPages;
using Microsoft.EntityFrameworkCore;
using System.Collections.Generic;
using System.Linq;

namespace EarlyBirdnBird.Pages.Host
{
    public class CalendarModel : PageModel
    {
        public string CurrentUserId { get; private set; }
        private readonly ApplicationDbContext _context;

        public CalendarModel(ApplicationDbContext context)
        {
            _context = context;
        }

        public List<ReservationWithRental> ReservationsWithRentals { get; private set; }

        public IActionResult OnGet()
        {
            CurrentUserId = User.FindFirst(System.Security.Claims.ClaimTypes.NameIdentifier)?.Value;

            ReservationsWithRentals = (from reservation in _context.Reservation
                                       join rental in _context.Rentals
                                       on Convert.ToInt32(reservation.Rental) equals rental.Id
                                       where rental.OwnerName.Id == CurrentUserId
                                       where reservation.Status != 3
                                       select new ReservationWithRental
                                       {
                                           Reservation = reservation,
                                           Rental = rental
                                       }).ToList();

            return Page();
        }
    }

    public class ReservationWithRental
    {
        public Reservation Reservation { get; set; }
        public Rental Rental { get; set; }
    }
}